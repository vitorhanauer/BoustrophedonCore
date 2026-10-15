#include "internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static int checked_array_size(size_t count, size_t element_size, size_t *bytes) {
    if (element_size != 0 && count > SIZE_MAX / element_size) return 0;
    *bytes = count * element_size;
    return 1;
}

static size_t disjoint_set_root(size_t *parents, size_t node) {
    if (parents[node] != node) parents[node] = disjoint_set_root(parents, parents[node]);
    return parents[node];
}

BpStatus bp_check_connected(const BpDecomposition *decomposition) {
    size_t *parents = NULL;
    unsigned char *ranks = NULL;
    size_t bytes;
    BpStatus status = BP_OK;

    if (!decomposition || decomposition->count == 0 || !decomposition->cells ||
        (decomposition->adjacency_count && !decomposition->adjacencies)) {
        return BP_ERR_ARGUMENT;
    }
    if (!checked_array_size(decomposition->count, sizeof(*parents), &bytes)) {
        return BP_ERR_LIMIT;
    }
    parents = malloc(bytes);
    if (!parents) return BP_ERR_MEMORY;
    ranks = calloc(decomposition->count, sizeof(*ranks));
    if (!ranks) {
        free(parents);
        return BP_ERR_MEMORY;
    }
    for (size_t i = 0; i < decomposition->count; ++i) parents[i] = i;

    for (size_t edge = 0; edge < decomposition->adjacency_count; ++edge) {
        size_t left = decomposition->adjacencies[edge].from;
        size_t right = decomposition->adjacencies[edge].to;
        size_t left_root, right_root;
        if (left >= decomposition->count || right >= decomposition->count) {
            status = BP_ERR_ARGUMENT;
            goto done;
        }
        left_root = disjoint_set_root(parents, left);
        right_root = disjoint_set_root(parents, right);
        if (left_root != right_root) {
            if (ranks[left_root] < ranks[right_root]) {
                parents[left_root] = right_root;
            } else if (ranks[left_root] > ranks[right_root]) {
                parents[right_root] = left_root;
            } else {
                parents[right_root] = left_root;
                ++ranks[left_root];
            }
        }
    }

    {
        size_t root = disjoint_set_root(parents, 0);
        for (size_t i = 1; i < decomposition->count; ++i) {
            size_t current = disjoint_set_root(parents, i);
            if (current != root) {
                status = BP_ERR_NO_PATH;
                break;
            }
        }
    }

done:
    free(parents);
    free(ranks);
    return status;
}

static BpStatus connector_distance(const BpDomain *domain, BpPoint from, BpPoint to,
        double *distance) {
    double result = 0.0;
    BpStatus status = bp_shortest_path(domain, from, to, NULL, &result);
    if (status != BP_OK) return status;
    if (!isfinite(result) || result < 0.0) return BP_ERR_LIMIT;
    *distance = result;
    return BP_OK;
}

static BpPoint path_first(const BpCell *cell, size_t variant) {
    return cell->variants[variant].points[0];
}

static BpPoint path_last(const BpCell *cell, size_t variant) {
    const BpPath *path = &cell->variants[variant];
    return path->points[path->count - 1];
}

/* Optimize orientations for a fixed cell order using a four-state DP. */
static BpStatus optimize_order(const BpDomain *domain, const BpDecomposition *decomposition,
        const size_t *order, BpPoint start, BpPoint end, const double *path_lengths,
        size_t *predecessors, size_t *variants, double *total_length) {
    double previous[4], current[4];
    size_t predecessor_index = 0;
    size_t first_cell = order[0];

    for (size_t variant = 0; variant < 4; ++variant) {
        double entry;
        BpStatus status = connector_distance(domain, start,
            path_first(&decomposition->cells[first_cell], variant), &entry);
        previous[variant] = status == BP_OK ? entry + path_lengths[first_cell * 4 + variant] : INFINITY;
        if (status != BP_OK && status != BP_ERR_NO_PATH) return status;
        if (!isfinite(previous[variant])) previous[variant] = INFINITY;
    }

    for (size_t position = 1; position < decomposition->count; ++position) {
        size_t previous_cell = order[position - 1];
        size_t cell = order[position];
        for (size_t variant = 0; variant < 4; ++variant) {
            double best = INFINITY;
            size_t best_previous = SIZE_MAX;
            for (size_t prior_variant = 0; prior_variant < 4; ++prior_variant) {
                double connector;
                BpStatus status;
                if (!isfinite(previous[prior_variant])) continue;
                status = connector_distance(domain,
                    path_last(&decomposition->cells[previous_cell], prior_variant),
                    path_first(&decomposition->cells[cell], variant), &connector);
                if (status == BP_ERR_NO_PATH) continue;
                if (status != BP_OK) return status;
                double candidate = previous[prior_variant] + connector + path_lengths[cell * 4 + variant];
                if (isfinite(candidate) && candidate < best) {
                    best = candidate;
                    best_previous = prior_variant;
                }
            }
            current[variant] = best;
            predecessors[predecessor_index++] = best_previous;
        }
        for (size_t variant = 0; variant < 4; ++variant) previous[variant] = current[variant];
    }

    {
        size_t last_cell = order[decomposition->count - 1];
        size_t best_variant = SIZE_MAX;
        double best = INFINITY;
        for (size_t variant = 0; variant < 4; ++variant) {
            double exit_cost;
            BpStatus status;
            if (!isfinite(previous[variant])) continue;
            status = connector_distance(domain, path_last(&decomposition->cells[last_cell], variant), end,
                &exit_cost);
            if (status == BP_ERR_NO_PATH) continue;
            if (status != BP_OK) return status;
            double candidate = previous[variant] + exit_cost;
            if (isfinite(candidate) && candidate < best) {
                best = candidate;
                best_variant = variant;
            }
        }
        if (best_variant == SIZE_MAX) return BP_ERR_NO_PATH;
        variants[decomposition->count - 1] = best_variant;
        for (size_t position = decomposition->count - 1; position > 0; --position) {
            size_t prior = predecessors[(position - 1) * 4 + variants[position]];
            if (prior == SIZE_MAX) return BP_ERR_NO_PATH;
            variants[position - 1] = prior;
        }
        *total_length = best;
    }
    return BP_OK;
}

static BpStatus greedy_order(const BpDomain *domain, const BpDecomposition *decomposition,
        BpPoint start, BpPoint end, const double *path_lengths, size_t *order) {
    size_t count = decomposition->count;
    unsigned char *visited = calloc(count, sizeof(*visited));
    if (!visited) return BP_ERR_MEMORY;

    BpPoint current = start;
    for (size_t position = 0; position < count; ++position) {
        size_t selected_cell = SIZE_MAX;
        size_t selected_variant = SIZE_MAX;
        double selected_cost = INFINITY;

        for (size_t cell = 0; cell < count; ++cell) {
            if (visited[cell]) continue;
            for (size_t variant = 0; variant < 4; ++variant) {
                double entry, exit_cost;
                BpStatus status = connector_distance(domain, current,
                    path_first(&decomposition->cells[cell], variant), &entry);
                if (status == BP_ERR_NO_PATH) continue;
                if (status != BP_OK) { free(visited); return status; }
                status = connector_distance(domain, path_last(&decomposition->cells[cell], variant), end,
                    &exit_cost);
                if (status == BP_ERR_NO_PATH) continue;
                if (status != BP_OK) { free(visited); return status; }
                double cost = entry + path_lengths[cell * 4 + variant] + exit_cost;
                if (!isfinite(cost)) continue;
                if (cost < selected_cost ||
                    (cost == selected_cost && (cell < selected_cell ||
                    (cell == selected_cell && variant < selected_variant)))) {
                    selected_cost = cost;
                    selected_cell = cell;
                    selected_variant = variant;
                }
            }
        }
        if (selected_cell == SIZE_MAX) { free(visited); return BP_ERR_NO_PATH; }
        order[position] = selected_cell;
        visited[selected_cell] = 1;
        current = path_last(&decomposition->cells[selected_cell], selected_variant);
    }
    free(visited);
    return BP_OK;
}

BpStatus bp_order_cells(const BpDomain *domain, const BpDecomposition *decomposition,
        BpPoint start, BpPoint end, BpVisitPlan *output) {
    size_t count, path_bytes, order_bytes, predecessor_count, predecessor_bytes;
    size_t *candidate_order = NULL, *predecessors = NULL, *candidate_variants = NULL;
    double *path_lengths = NULL;
    BpStatus status;

    if (!output) return BP_ERR_ARGUMENT;
    *output = (BpVisitPlan){0};
    if (!domain || !decomposition || decomposition->count == 0 || !decomposition->cells ||
        !isfinite(start.x) || !isfinite(start.y) || !isfinite(end.x) || !isfinite(end.y)) {
        return BP_ERR_ARGUMENT;
    }
    count = decomposition->count;
    status = bp_check_connected(decomposition);
    if (status != BP_OK) return status;
    if (!checked_array_size(count, sizeof(*output->order), &order_bytes) ||
        !checked_array_size(count, 4, &predecessor_count) ||
        !checked_array_size(predecessor_count, sizeof(*predecessors), &predecessor_bytes) ||
        !checked_array_size(predecessor_count, sizeof(*path_lengths), &path_bytes)) {
        return BP_ERR_LIMIT;
    }
    output->order = malloc(order_bytes);
    output->variants = malloc(order_bytes);
    path_lengths = malloc(path_bytes);
    predecessors = malloc(predecessor_bytes);
    candidate_order = malloc(order_bytes);
    candidate_variants = malloc(order_bytes);
    if (!output->order || !output->variants || !path_lengths || !predecessors ||
        !candidate_order || !candidate_variants) {
        status = BP_ERR_MEMORY;
        goto cleanup;
    }
    output->count = count;

    for (size_t cell = 0; cell < count; ++cell) {
        for (size_t variant = 0; variant < 4; ++variant) {
            const BpPath *path = &decomposition->cells[cell].variants[variant];
            if (!path->points || path->count == 0) { status = BP_ERR_ARGUMENT; goto cleanup; }
            double length = bp_path_length(path);
            if (!isfinite(length) || length < 0.0) { status = BP_ERR_LIMIT; goto cleanup; }
            path_lengths[cell * 4 + variant] = length;
        }
    }

    status = greedy_order(domain, decomposition, start, end, path_lengths, output->order);
    if (status != BP_OK) goto cleanup;
    status = optimize_order(domain, decomposition, output->order, start, end,
        path_lengths, predecessors, output->variants, &output->initial_length);
    if (status != BP_OK) goto cleanup;
    output->length = output->initial_length;

    for (size_t pass = 0; pass < 10; ++pass) {
        int improved = 0;
        double best_length = output->length;
        for (size_t first = 0; first + 1 < count; ++first) {
            for (size_t last = first + 1; last < count; ++last) {
                for (size_t i = 0; i < count; ++i) candidate_order[i] = output->order[i];
                for (size_t left = first, right = last; left < right; ++left, --right) {
                    size_t temporary = candidate_order[left];
                    candidate_order[left] = candidate_order[right];
                    candidate_order[right] = temporary;
                }
                double candidate_length;
                status = optimize_order(domain, decomposition, candidate_order, start, end,
                    path_lengths, predecessors, candidate_variants, &candidate_length);
                if (status == BP_ERR_NO_PATH) continue;
                if (status != BP_OK) goto cleanup;
                if (candidate_length + 1e-12 < best_length) {
                    best_length = candidate_length;
                    for (size_t i = 0; i < count; ++i) {
                        output->order[i] = candidate_order[i];
                        output->variants[i] = candidate_variants[i];
                    }
                    improved = 1;
                }
            }
        }
        output->length = best_length;
        if (!improved) break;
    }
    status = BP_OK;

cleanup:
    free(path_lengths);
    free(predecessors);
    free(candidate_order);
    free(candidate_variants);
    if (status != BP_OK) bp_free_plan(output);
    return status;
}
