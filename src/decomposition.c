#include "internal.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct { BpPoint a, b; } SweepEdge;
typedef struct {
    long double y;
    SweepEdge edge;
} Crossing;
typedef struct {
    SweepEdge lower_edge, upper_edge;
    BpPoint lower_left, lower_right, upper_left, upper_right;
    size_t cell;
} SlabInterval;
typedef struct {
    SlabInterval *items;
    size_t count, capacity;
} IntervalList;
typedef struct {
    BpPoint *points;
    size_t count, capacity;
} PointList;
typedef struct {
    PointList lower, upper;
} CellBuild;
typedef struct {
    BpAdjacency *items;
    size_t count, capacity;
} AdjacencyList;

static int compare_double(const void *left, const void *right) {
    double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

static int compare_crossing(const void *left, const void *right) {
    long double a = ((const Crossing *)left)->y;
    long double b = ((const Crossing *)right)->y;
    return (a > b) - (a < b);
}

static BpStatus reserve(void **memory, size_t *capacity, size_t needed, size_t item_size) {
    size_t next;
    void *replacement;
    if (needed <= *capacity) return BP_OK;
    if (needed > SIZE_MAX / item_size) return BP_ERR_LIMIT;
    next = *capacity ? *capacity : 8;
    while (next < needed) {
        if (next > SIZE_MAX / 2) { next = needed; break; }
        next *= 2;
    }
    if (next > SIZE_MAX / item_size) return BP_ERR_LIMIT;
    replacement = realloc(*memory, next * item_size);
    if (!replacement) return BP_ERR_MEMORY;
    *memory = replacement;
    *capacity = next;
    return BP_OK;
}

static BpStatus append_point(PointList *list, BpPoint point) {
    BpStatus status;
    if (list->count && bp_same_point(list->points[list->count - 1], point)) return BP_OK;
    status = reserve((void **)&list->points, &list->capacity, list->count + 1, sizeof(*list->points));
    if (status != BP_OK) return status;
    list->points[list->count++] = point;
    return BP_OK;
}

static long double edge_y_at(const SweepEdge *edge, double x) {
    long double dx = (long double)edge->b.x - edge->a.x;
    if (dx == 0) return (long double)edge->a.y;
    return (long double)edge->a.y +
        (((long double)x - edge->a.x) / dx) * ((long double)edge->b.y - edge->a.y);
}

static BpPoint point_on_edge(const SweepEdge *edge, double x) {
    BpPoint point = { x, (double)edge_y_at(edge, x) };
    return point;
}

static int point_is_free(const BpDomain *domain, BpPoint point) {
    if (bp_point_location(domain->area, point) < 0) return 0;
    for (size_t i = 0; i < domain->obstacle_count; ++i)
        if (bp_point_location(domain->obstacles[i], point) > 0) return 0;
    return 1;
}

static BpStatus add_ring_crossings(BpPolygon ring, double x, Crossing *crossings, size_t *count) {
    size_t n = bp_ring_count(ring);
    for (size_t i = 0; i < n; ++i) {
        BpPoint a = ring.vertices[i], b = ring.vertices[(i + 1) % n];
        if ((a.x < x && x < b.x) || (b.x < x && x < a.x)) {
            SweepEdge edge = { a, b };
            crossings[*count].y = edge_y_at(&edge, x);
            crossings[*count].edge = edge;
            ++*count;
        }
    }
    return BP_OK;
}

static BpStatus slab_intervals(const BpDomain *domain, double left, double right,
    Crossing *crossings, size_t crossing_capacity, IntervalList *output) {
    double middle = left / 2.0 + right / 2.0;
    size_t count = 0;
    BpStatus status;
    status = add_ring_crossings(domain->area, middle, crossings, &count);
    if (status != BP_OK) return status;
    for (size_t i = 0; i < domain->obstacle_count; ++i) {
        status = add_ring_crossings(domain->obstacles[i], middle, crossings, &count);
        if (status != BP_OK) return status;
    }
    if (count > crossing_capacity) return BP_ERR_LIMIT;
    qsort(crossings, count, sizeof(*crossings), compare_crossing);
    for (size_t i = 0; i + 1 < count; ++i) {
        long double low = crossings[i].y, high = crossings[i + 1].y;
        if (!(low < high)) continue;
        BpPoint sample = { middle, (double)(low / 2.0L + high / 2.0L) };
        if (point_is_free(domain, sample)) {
            SlabInterval interval;
            interval.lower_edge = crossings[i].edge;
            interval.upper_edge = crossings[i + 1].edge;
            interval.lower_left = point_on_edge(&interval.lower_edge, left);
            interval.lower_right = point_on_edge(&interval.lower_edge, right);
            interval.upper_left = point_on_edge(&interval.upper_edge, left);
            interval.upper_right = point_on_edge(&interval.upper_edge, right);
            interval.cell = SIZE_MAX;
            status = reserve((void **)&output->items, &output->capacity,
                output->count + 1, sizeof(*output->items));
            if (status != BP_OK) return status;
            output->items[output->count++] = interval;
        }
    }
    return BP_OK;
}

static int positive_overlap_at(const SlabInterval *left, const SlabInterval *right, double x) {
    long double left_low = edge_y_at(&left->lower_edge, x);
    long double left_high = edge_y_at(&left->upper_edge, x);
    long double right_low = edge_y_at(&right->lower_edge, x);
    long double right_high = edge_y_at(&right->upper_edge, x);
    long double low = left_low > right_low ? left_low : right_low;
    long double high = left_high < right_high ? left_high : right_high;
    return high > low;
}

static BpStatus append_cell(CellBuild **cells, size_t *count, size_t *capacity,
    SlabInterval *interval, size_t *cell_index) {
    BpStatus status = reserve((void **)cells, capacity, *count + 1, sizeof(**cells));
    if (status != BP_OK) return status;
    (*cells)[*count] = (CellBuild){0};
    *cell_index = *count;
    ++*count;
    interval->cell = *cell_index;
    status = append_point(&(*cells)[*cell_index].lower, interval->lower_left);
    if (status == BP_OK) status = append_point(&(*cells)[*cell_index].upper, interval->upper_left);
    if (status == BP_OK) status = append_point(&(*cells)[*cell_index].lower, interval->lower_right);
    if (status == BP_OK) status = append_point(&(*cells)[*cell_index].upper, interval->upper_right);
    return status;
}

static BpStatus extend_cell(CellBuild *cell, SlabInterval *interval) {
    BpStatus status = append_point(&cell->lower, interval->lower_left);
    if (status == BP_OK) status = append_point(&cell->lower, interval->lower_right);
    if (status == BP_OK) status = append_point(&cell->upper, interval->upper_left);
    if (status == BP_OK) status = append_point(&cell->upper, interval->upper_right);
    return status;
}

static BpStatus add_adjacency(AdjacencyList *list, size_t a, size_t b) {
    BpStatus status;
    if (a == b) return BP_OK;
    if (a > b) { size_t swap = a; a = b; b = swap; }
    for (size_t i = 0; i < list->count; ++i)
        if (list->items[i].from == a && list->items[i].to == b) return BP_OK;
    status = reserve((void **)&list->items, &list->capacity,
        list->count + 1, sizeof(*list->items));
    if (status != BP_OK) return status;
    list->items[list->count++] = (BpAdjacency){ a, b };
    return BP_OK;
}

static BpStatus finish_cell(BpCell *cell, const CellBuild *build) {
    BpStatus status;
    if (build->lower.count + build->upper.count > SIZE_MAX / sizeof(BpPoint)) return BP_ERR_LIMIT;
    if (build->lower.count < 2 || build->upper.count < 2) return BP_ERR_GEOMETRY;
    for (size_t i = 0; i < build->lower.count; ++i) {
        status = bp_path_append(&cell->boundary, build->lower.points[i]);
        if (status != BP_OK) return status;
    }
    for (size_t i = build->upper.count; i > 0; --i) {
        status = bp_path_append(&cell->boundary, build->upper.points[i - 1]);
        if (status != BP_OK) return status;
    }
    if (cell->boundary.count > 1 && bp_same_point(cell->boundary.points[0],
        cell->boundary.points[cell->boundary.count - 1])) --cell->boundary.count;
    return cell->boundary.count >= 3 ? BP_OK : BP_ERR_GEOMETRY;
}

static void free_builds(CellBuild *cells, size_t count) {
    if (!cells) return;
    for (size_t i = 0; i < count; ++i) {
        free(cells[i].lower.points);
        free(cells[i].upper.points);
    }
    free(cells);
}

BpStatus bp_decompose(const BpDomain *domain, BpDecomposition *output) {
    double *events = NULL;
    Crossing *crossings = NULL;
    IntervalList previous = {0}, current = {0};
    CellBuild *builds = NULL;
    size_t event_count = 0, event_capacity = 0, crossing_capacity = 0;
    size_t cell_count = 0, cell_capacity = 0;
    AdjacencyList adjacencies = {0};
    BpDecomposition result = {0};
    BpStatus status = BP_OK;

    if (!domain || !output || !domain->area.vertices || domain->area.count < 3 ||
        (domain->obstacle_count && !domain->obstacles)) return BP_ERR_ARGUMENT;
    *output = (BpDecomposition){0};

    size_t total_vertices = bp_ring_count(domain->area);
    for (size_t i = 0; i < domain->obstacle_count; ++i) {
        size_t n = bp_ring_count(domain->obstacles[i]);
        if (!domain->obstacles[i].vertices || n < 3 || n > SIZE_MAX - total_vertices)
            return BP_ERR_ARGUMENT;
        total_vertices += n;
    }
    if (total_vertices < 3 || total_vertices > SIZE_MAX / sizeof(Crossing) ||
        total_vertices > SIZE_MAX / sizeof(double)) return BP_ERR_LIMIT;
    crossing_capacity = total_vertices;
    crossings = malloc(crossing_capacity * sizeof(*crossings));
    if (!crossings) return BP_ERR_MEMORY;

    for (size_t ring_index = 0; ring_index <= domain->obstacle_count; ++ring_index) {
        BpPolygon ring = ring_index ? domain->obstacles[ring_index - 1] : domain->area;
        size_t n = bp_ring_count(ring);
        for (size_t i = 0; i < n; ++i) {
            if (!isfinite(ring.vertices[i].x) || !isfinite(ring.vertices[i].y)) {
                status = BP_ERR_ARGUMENT; goto cleanup;
            }
            status = reserve((void **)&events, &event_capacity, event_count + 1, sizeof(*events));
            if (status != BP_OK) goto cleanup;
            events[event_count++] = ring.vertices[i].x;
        }
    }
    qsort(events, event_count, sizeof(*events), compare_double);
    size_t unique_events = 0;
    for (size_t i = 0; i < event_count; ++i)
        if (!unique_events || events[i] != events[unique_events - 1]) events[unique_events++] = events[i];
    event_count = unique_events;
    if (event_count < 2) { status = BP_ERR_GEOMETRY; goto cleanup; }

    for (size_t slab = 0; slab + 1 < event_count; ++slab) {
        double left = events[slab], right = events[slab + 1];
        if (!(left < right)) continue;
        current.count = 0;
        status = slab_intervals(domain, left, right, crossings, crossing_capacity, &current);
        if (status != BP_OK) goto cleanup;

        if (!previous.count) {
            for (size_t j = 0; j < current.count; ++j) {
                size_t ignored;
                status = append_cell(&builds, &cell_count, &cell_capacity, &current.items[j], &ignored);
                if (status != BP_OK) goto cleanup;
            }
        } else {
            size_t *previous_degrees = calloc(previous.count, sizeof(*previous_degrees));
            size_t *current_degrees = calloc(current.count, sizeof(*current_degrees));
            if (!previous_degrees || !current_degrees) {
                free(previous_degrees); free(current_degrees); status = BP_ERR_MEMORY; goto cleanup;
            }
            for (size_t i = 0; i < previous.count; ++i) {
                for (size_t j = 0; j < current.count; ++j) {
                    if (positive_overlap_at(&previous.items[i], &current.items[j], left)) {
                        ++previous_degrees[i]; ++current_degrees[j];
                    }
                }
            }
            for (size_t j = 0; j < current.count; ++j) {
                size_t predecessor = SIZE_MAX, predecessor_count = 0;
                for (size_t i = 0; i < previous.count; ++i) {
                    if (positive_overlap_at(&previous.items[i], &current.items[j], left)) {
                        predecessor = i;
                        ++predecessor_count;
                    }
                }
                if (predecessor_count == 1 && current_degrees[j] == 1 &&
                    previous_degrees[predecessor] == 1) {
                    current.items[j].cell = previous.items[predecessor].cell;
                    status = extend_cell(&builds[current.items[j].cell], &current.items[j]);
                } else {
                    size_t created;
                    status = append_cell(&builds, &cell_count, &cell_capacity, &current.items[j], &created);
                }
                if (status != BP_OK) { free(previous_degrees); free(current_degrees); goto cleanup; }
            }
            for (size_t i = 0; i < previous.count; ++i) {
                for (size_t j = 0; j < current.count; ++j) {
                    if (positive_overlap_at(&previous.items[i], &current.items[j], left) &&
                        previous.items[i].cell != current.items[j].cell) {
                        status = add_adjacency(&adjacencies, previous.items[i].cell, current.items[j].cell);
                        if (status != BP_OK) { free(previous_degrees); free(current_degrees); goto cleanup; }
                    }
                }
            }
            free(previous_degrees);
            free(current_degrees);
        }
        SlabInterval *swap_items = previous.items;
        size_t swap_capacity = previous.capacity;
        previous.items = current.items; previous.count = current.count; previous.capacity = current.capacity;
        current.items = swap_items; current.count = 0; current.capacity = swap_capacity;
    }

    if (!cell_count || cell_count > SIZE_MAX / sizeof(*result.cells) ||
        adjacencies.count > SIZE_MAX / sizeof(*result.adjacencies)) {
        status = BP_ERR_LIMIT; goto cleanup;
    }
    result.cells = calloc(cell_count, sizeof(*result.cells));
    if (!result.cells) { status = BP_ERR_MEMORY; goto cleanup; }
    result.count = cell_count;
    for (size_t i = 0; i < cell_count; ++i) {
        status = finish_cell(&result.cells[i], &builds[i]);
        if (status != BP_OK) { bp_free_decomposition(&result); goto cleanup; }
    }
    if (adjacencies.count) {
        result.adjacencies = adjacencies.items;
        result.adjacency_count = adjacencies.count;
        adjacencies.items = NULL;
    }
    *output = result;
    result = (BpDecomposition){0};

cleanup:
    bp_free_decomposition(&result);
    free(events); free(crossings);
    free(previous.items); free(current.items);
    free(adjacencies.items);
    free_builds(builds, cell_count);
    return status;
}
