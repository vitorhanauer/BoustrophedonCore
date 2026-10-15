#include "internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static BpStatus append_path(BpPath *destination, const BpPath *source) {
    if (!source || (source->count && !source->points)) return BP_ERR_ARGUMENT;
    for (size_t i = 0; i < source->count; ++i) {
        BpStatus status = bp_path_append(destination, source->points[i]);
        if (status != BP_OK) return status;
    }
    return BP_OK;
}

static BpStatus shortest_distance(BpPoint a, BpPoint b, double *result) {
    double length = hypot(b.x - a.x, b.y - a.y);
    if (!isfinite(length)) return BP_ERR_LIMIT;
    *result = length;
    return BP_OK;
}

static BpStatus point_index(const BpPoint *points, size_t count, BpPoint point,
        size_t *index) {
    for (size_t i = 0; i < count; ++i) {
        if (bp_same_point(points[i], point)) {
            *index = i;
            return BP_OK;
        }
    }
    if (count >= SIZE_MAX / sizeof(*points)) return BP_ERR_LIMIT;
    *index = count;
    return BP_OK;
}

/* Dijkstra over the implicit visibility graph.  Only the current distance,
 * parent, and visited state are stored, so memory use stays linear in vertices. */
BpStatus bp_shortest_path(const BpDomain *domain, BpPoint start, BpPoint end,
        BpPath *output, double *distance) {
    if (distance) *distance = 0.0;
    if (!domain || !domain->area.vertices || domain->area.count < 3
        || (domain->obstacle_count && !domain->obstacles)
        || !isfinite(start.x) || !isfinite(start.y)
        || !isfinite(end.x) || !isfinite(end.y)) return BP_ERR_ARGUMENT;
    if (output && (output->points || output->count)) return BP_ERR_ARGUMENT;
    if (!bp_point_free(domain, start) || !bp_point_free(domain, end)) return BP_ERR_GEOMETRY;

    int visible = 0;
    BpStatus status = bp_segment_visible(domain, start, end, &visible);
    if (status != BP_OK) return status;
    if (visible) {
        double length;
        status = shortest_distance(start, end, &length);
        if (status != BP_OK) return status;
        if (output) {
            BpPath direct = {0};
            status = bp_path_append(&direct, start);
            if (status == BP_OK) status = bp_path_append(&direct, end);
            if (status != BP_OK) { freeWaypoints(&direct); return status; }
            *output = direct;
        }
        if (distance) *distance = length;
        return BP_OK;
    }

    size_t vertex_count = bp_ring_count(domain->area);
    for (size_t i = 0; i < domain->obstacle_count; ++i) {
        size_t ring_count = bp_ring_count(domain->obstacles[i]);
        if (!domain->obstacles[i].vertices || ring_count < 3
            || ring_count > SIZE_MAX - vertex_count) return BP_ERR_ARGUMENT;
        vertex_count += ring_count;
    }
    if (vertex_count > SIZE_MAX - 2 || vertex_count + 2 > SIZE_MAX / sizeof(BpPoint))
        return BP_ERR_LIMIT;
    size_t capacity = vertex_count + 2;
    BpPoint *vertices = malloc(capacity * sizeof(*vertices));
    if (!vertices) return BP_ERR_MEMORY;
    size_t count = 0, start_index = SIZE_MAX, end_index = SIZE_MAX;
    vertices[count++] = start;
    if (!bp_same_point(start, end)) {
        vertices[count++] = end;
        end_index = 1;
    } else {
        end_index = 0;
    }
    for (size_t ring_index = 0; ring_index <= domain->obstacle_count; ++ring_index) {
        BpPolygon ring = ring_index ? domain->obstacles[ring_index - 1] : domain->area;
        size_t n = bp_ring_count(ring);
        for (size_t i = 0; i < n; ++i) {
            size_t index;
            status = point_index(vertices, count, ring.vertices[i], &index);
            if (status != BP_OK) { free(vertices); return status; }
            if (index == count) vertices[count++] = ring.vertices[i];
        }
    }
    start_index = 0;
    if (bp_same_point(start, end)) end_index = 0;

    if (count > SIZE_MAX / sizeof(double) || count > SIZE_MAX / sizeof(size_t)) {
        free(vertices);
        return BP_ERR_LIMIT;
    }
    double *best = malloc(count * sizeof(*best));
    size_t *parent = malloc(count * sizeof(*parent));
    unsigned char *settled = calloc(count, sizeof(*settled));
    if (!best || !parent || !settled) {
        free(vertices); free(best); free(parent); free(settled);
        return BP_ERR_MEMORY;
    }
    for (size_t i = 0; i < count; ++i) {
        best[i] = INFINITY;
        parent[i] = SIZE_MAX;
    }
    best[start_index] = 0.0;

    for (size_t iteration = 0; iteration < count; ++iteration) {
        size_t current = SIZE_MAX;
        for (size_t i = 0; i < count; ++i) {
            if (!settled[i] && (current == SIZE_MAX || best[i] < best[current])) current = i;
        }
        if (current == SIZE_MAX || !isfinite(best[current])) break;
        if (current == end_index) break;
        settled[current] = 1;
        for (size_t next = 0; next < count; ++next) {
            if (settled[next] || next == current) continue;
            status = bp_segment_visible(domain, vertices[current], vertices[next], &visible);
            if (status != BP_OK) goto cleanup;
            if (!visible) continue;
            double edge_length;
            status = shortest_distance(vertices[current], vertices[next], &edge_length);
            if (status != BP_OK) goto cleanup;
            double candidate = best[current] + edge_length;
            if (!isfinite(candidate)) { status = BP_ERR_LIMIT; goto cleanup; }
            if (candidate < best[next]) {
                best[next] = candidate;
                parent[next] = current;
            }
        }
    }
    if (!isfinite(best[end_index])) { status = BP_ERR_NO_PATH; goto cleanup; }

    if (output) {
        size_t path_count = 1;
        for (size_t index = end_index; index != start_index; index = parent[index]) {
            if (parent[index] == SIZE_MAX || path_count >= count) {
                status = BP_ERR_NO_PATH;
                goto cleanup;
            }
            ++path_count;
        }
        size_t *indices = malloc(path_count * sizeof(*indices));
        if (!indices) { status = BP_ERR_MEMORY; goto cleanup; }
        size_t index = end_index;
        for (size_t position = path_count; position-- > 0;) {
            indices[position] = index;
            if (index == start_index) break;
            index = parent[index];
        }
        BpPath route = {0};
        for (size_t i = 0; i < path_count; ++i) {
            status = bp_path_append(&route, vertices[indices[i]]);
            if (status != BP_OK) break;
        }
        free(indices);
        if (status != BP_OK) { freeWaypoints(&route); goto cleanup; }
        *output = route;
    }
    if (distance) *distance = best[end_index];
    status = BP_OK;

cleanup:
    free(vertices); free(best); free(parent); free(settled);
    return status;
}

BpStatus bp_stitch(const BpDomain *domain, const BpDecomposition *decomposition,
        const BpVisitPlan *plan, BpPoint start, BpPoint end, BpPath *output) {
    if (!domain || !decomposition || !plan || !output || !domain->area.vertices
        || (domain->obstacle_count && !domain->obstacles)
        || !isfinite(start.x) || !isfinite(start.y)
        || !isfinite(end.x) || !isfinite(end.y)) return BP_ERR_ARGUMENT;
    if (output->points || output->count) return BP_ERR_ARGUMENT;
    if (plan->count != decomposition->count
        || (plan->count && (!plan->order || !plan->variants || !decomposition->cells)))
        return BP_ERR_ARGUMENT;
    if (!bp_point_free(domain, start) || !bp_point_free(domain, end)) return BP_ERR_GEOMETRY;
    if (plan->count > SIZE_MAX / sizeof(unsigned char)) return BP_ERR_LIMIT;
    unsigned char *seen = decomposition->count ? calloc(decomposition->count, 1) : NULL;
    if (decomposition->count && !seen) return BP_ERR_MEMORY;
    for (size_t i = 0; i < plan->count; ++i) {
        size_t cell_index = plan->order[i], variant_index = plan->variants[i];
        if (cell_index >= decomposition->count || variant_index >= 4 || seen[cell_index]) {
            free(seen);
            return BP_ERR_ARGUMENT;
        }
        seen[cell_index] = 1;
        const BpPath *path = &decomposition->cells[cell_index].variants[variant_index];
        if (!path->points || path->count < 2) { free(seen); return BP_ERR_GEOMETRY; }
        for (size_t point = 0; point < path->count; ++point) {
            if (!isfinite(path->points[point].x) || !isfinite(path->points[point].y)
                || !bp_point_free(domain, path->points[point])) {
                free(seen);
                return BP_ERR_GEOMETRY;
            }
            if (point) {
                int visible = 0;
                BpStatus status = bp_segment_visible(domain, path->points[point - 1],
                    path->points[point], &visible);
                if (status != BP_OK) { free(seen); return status; }
                if (!visible) { free(seen); return BP_ERR_GEOMETRY; }
            }
        }
    }
    free(seen);

    BpPath route = {0};
    BpStatus status = bp_path_append(&route, start);
    BpPoint current = start;
    for (size_t i = 0; status == BP_OK && i < plan->count; ++i) {
        const BpPath *cell_path = &decomposition->cells[plan->order[i]].variants[plan->variants[i]];
        BpPath connector = {0};
        status = bp_shortest_path(domain, current, cell_path->points[0], &connector, NULL);
        if (status == BP_OK) status = append_path(&route, &connector);
        freeWaypoints(&connector);
        if (status == BP_OK) status = append_path(&route, cell_path);
        if (status == BP_OK) current = cell_path->points[cell_path->count - 1];
    }
    if (status == BP_OK) {
        BpPath connector = {0};
        status = bp_shortest_path(domain, current, end, &connector, NULL);
        if (status == BP_OK) status = append_path(&route, &connector);
        freeWaypoints(&connector);
    }
    if (status == BP_OK && !route.count) status = BP_ERR_NO_PATH;
    if (status == BP_OK && !isfinite(bp_path_length(&route))) status = BP_ERR_LIMIT;
    if (status != BP_OK) { freeWaypoints(&route); return status; }
    *output = route;
    return BP_OK;
}
