#include "test_support.h"

static const BpPoint rectangle[] = {
    {0, 0}, {10, 0}, {10, 10}, {0, 10}
};
static const BpPoint triangle[] = {
    {0, 0}, {10, 0}, {5, 10}
};
static const BpPoint concave[] = {
    {0, 0}, {4, 0}, {4, 1}, {1, 1},
    {1, 3}, {4, 3}, {4, 4}, {0, 4}
};
static const BpPoint obstacle_vertices[] = {
    {4, 3}, {6, 3}, {6, 7}, {4, 7}
};

static double ring_area(const BpPoint *vertices, size_t count) {
    size_t usable = count;
    double twice_area = 0.0;
    size_t i;
    if (usable > 1 && bp_same_point(vertices[0], vertices[usable - 1])) {
        --usable;
    }
    for (i = 0; i < usable; ++i) {
        const BpPoint a = vertices[i];
        const BpPoint b = vertices[(i + 1) % usable];
        twice_area += a.x * b.y - b.x * a.y;
    }
    return fabs(twice_area) / 2.0;
}

static double cell_area(const BpCell *cell) {
    CHECK(cell->boundary.count >= 3);
    return ring_area(cell->boundary.points, cell->boundary.count);
}

static int edges_share_positive_length(const BpPath *a, const BpPath *b) {
    size_t i, j;
    for (i = 0; i < a->count; ++i) {
        BpPoint a0 = a->points[i];
        BpPoint a1 = a->points[(i + 1) % a->count];
        for (j = 0; j < b->count; ++j) {
            BpPoint b0 = b->points[j];
            BpPoint b1 = b->points[(j + 1) % b->count];
            double dx = a1.x - a0.x;
            double dy = a1.y - a0.y;
            double cross0 = dx * (b0.y - a0.y) - dy * (b0.x - a0.x);
            double cross1 = dx * (b1.y - a0.y) - dy * (b1.x - a0.x);
            double scale = fabs(dx) + fabs(dy) + 1.0;
            double overlap;
            if (fabs(cross0) > 1e-9 * scale || fabs(cross1) > 1e-9 * scale) {
                continue;
            }
            if (fabs(dx) >= fabs(dy)) {
                double low = fmax(fmin(a0.x, a1.x), fmin(b0.x, b1.x));
                double high = fmin(fmax(a0.x, a1.x), fmax(b0.x, b1.x));
                overlap = high - low;
            } else {
                double low = fmax(fmin(a0.y, a1.y), fmin(b0.y, b1.y));
                double high = fmin(fmax(a0.y, a1.y), fmax(b0.y, b1.y));
                overlap = high - low;
            }
            if (overlap > 1e-9) {
                return 1;
            }
        }
    }
    return 0;
}

static void check_adjacencies(const BpDecomposition *decomposition) {
    size_t i, j;
    size_t *queue;
    unsigned char *visited;
    if (decomposition->count <= 1) {
        CHECK(decomposition->adjacency_count == 0);
        return;
    }
    CHECK(decomposition->adjacency_count > 0);
    queue = calloc(decomposition->count, sizeof(*queue));
    visited = calloc(decomposition->count, sizeof(*visited));
    CHECK(queue != NULL && visited != NULL);
    for (i = 0; i < decomposition->adjacency_count; ++i) {
        const BpAdjacency edge = decomposition->adjacencies[i];
        CHECK(edge.from < decomposition->count);
        CHECK(edge.to < decomposition->count);
        CHECK(edge.from != edge.to);
        CHECK(edges_share_positive_length(&decomposition->cells[edge.from].boundary,
            &decomposition->cells[edge.to].boundary));
        for (j = 0; j < i; ++j) {
            const BpAdjacency previous = decomposition->adjacencies[j];
            CHECK(!((edge.from == previous.from && edge.to == previous.to) ||
                    (edge.from == previous.to && edge.to == previous.from)));
        }
    }
    {
        size_t head = 0, tail = 0;
        queue[tail++] = 0;
        visited[0] = 1;
        while (head < tail) {
            size_t current = queue[head++];
            for (i = 0; i < decomposition->adjacency_count; ++i) {
                BpAdjacency edge = decomposition->adjacencies[i];
                size_t next = decomposition->count;
                if (edge.from == current) next = edge.to;
                if (edge.to == current) next = edge.from;
                if (next < decomposition->count && !visited[next]) {
                    visited[next] = 1;
                    queue[tail++] = next;
                }
            }
        }
    }
    for (i = 0; i < decomposition->count; ++i) CHECK(visited[i]);
    free(queue);
    free(visited);
}

static BpDecomposition decompose(const BpPoint *area_vertices, size_t area_count,
    const BpPolygon *obstacles, size_t obstacle_count, size_t expected_cells,
    double expected_area) {
    BpDomain domain = {
        {area_vertices, area_count}, obstacles, obstacle_count
    };
    BpDecomposition result = {0};
    double total_area = 0.0;
    size_t i;
    CHECK(bp_decompose(&domain, &result) == BP_OK);
    CHECK(result.count == expected_cells);
    CHECK(result.cells != NULL);
    for (i = 0; i < result.count; ++i) {
        double area = cell_area(&result.cells[i]);
        CHECK(area > 1e-9);
        total_area += area;
    }
    NEAR(total_area, expected_area, 1e-8 * fmax(1.0, expected_area));
    check_adjacencies(&result);
    return result;
}

static void check_free(BpDecomposition *decomposition) {
    bp_free_decomposition(decomposition);
    CHECK(decomposition->cells == NULL);
    CHECK(decomposition->adjacencies == NULL);
    CHECK(decomposition->count == 0);
    CHECK(decomposition->adjacency_count == 0);
}

static void test_basic_shapes(void) {
    BpPolygon obstacle = {obstacle_vertices, sizeof(obstacle_vertices) / sizeof(obstacle_vertices[0])};
    BpDecomposition result;
    result = decompose(rectangle, 4, NULL, 0, 1, 100.0);
    check_free(&result);
    result = decompose(triangle, 3, NULL, 0, 1, 50.0);
    check_free(&result);
    result = decompose(concave, 8, NULL, 0, 3, 10.0);
    check_free(&result);
    result = decompose(rectangle, 4, &obstacle, 1, 4, 92.0);
    check_free(&result);
}

static void test_coincident_collinear_events(void) {
    /* Extra collinear vertices exercise event deduplication at x=4 and x=6. */
    const BpPoint area[] = {
        {0, 0}, {4, 0}, {6, 0}, {10, 0}, {10, 10},
        {6, 10}, {4, 10}, {0, 10}, {0, 0}
    };
    const BpPoint obstacle_ring[] = {
        {4, 3}, {5, 3}, {6, 3}, {6, 5}, {6, 7},
        {5, 7}, {4, 7}, {4, 5}, {4, 3}
    };
    BpPolygon obstacle = {obstacle_ring, sizeof(obstacle_ring) / sizeof(obstacle_ring[0])};
    BpDecomposition result = decompose(area,
        sizeof(area) / sizeof(area[0]), &obstacle, 1, 4, 92.0);
    check_free(&result);
}

static void test_orientation_closure_and_scale(void) {
    const BpPoint reverse_closed[] = {
        {0, 4}, {4, 4}, {4, 3}, {1, 3}, {1, 1},
        {4, 1}, {4, 0}, {0, 0}, {0, 4}
    };
    const BpPoint scaled[] = {
        {0, 0}, {40, 0}, {40, 10}, {10, 10},
        {10, 30}, {40, 30}, {40, 40}, {0, 40}
    };
    BpDecomposition result;
    result = decompose(reverse_closed,
        sizeof(reverse_closed) / sizeof(reverse_closed[0]), NULL, 0, 3, 10.0);
    check_free(&result);
    result = decompose(scaled, sizeof(scaled) / sizeof(scaled[0]), NULL, 0, 3, 1000.0);
    check_free(&result);
}

int main(void) {
    test_basic_shapes();
    test_coincident_collinear_events();
    test_orientation_closure_and_scale();
    return EXIT_SUCCESS;
}
