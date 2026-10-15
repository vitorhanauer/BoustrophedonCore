#include "test_support.h"

static BpPoint rectangle_vertices[] = {
    {0.0, 0.0}, {100.0, 0.0}, {100.0, 100.0}, {0.0, 100.0}
};

static BpDomain rectangle_domain(void) {
    BpDomain domain = {
        .area = {rectangle_vertices, sizeof(rectangle_vertices) / sizeof(rectangle_vertices[0])},
        .obstacles = NULL,
        .obstacle_count = 0
    };
    return domain;
}

static void set_cell_paths(BpCell *cell, BpPoint points[4]) {
    for (size_t variant = 0; variant < 4; ++variant) {
        cell->variants[variant].points = &points[variant];
        cell->variants[variant].count = 1;
    }
}

static void test_connectivity(void) {
    BpCell cells[4] = {0};
    BpAdjacency chain_edges[] = {{0, 1}, {1, 2}, {2, 3}};
    BpDecomposition chain = {cells, 4, chain_edges, 3};
    CHECK(bp_check_connected(&chain) == BP_OK);

    BpAdjacency branch_edges[] = {{0, 1}, {0, 2}, {2, 3}};
    BpDecomposition branch = {cells, 4, branch_edges, 3};
    CHECK(bp_check_connected(&branch) == BP_OK);

    BpAdjacency disconnected_edges[] = {{0, 1}, {2, 3}};
    BpDecomposition disconnected = {cells, 4, disconnected_edges, 2};
    CHECK(bp_check_connected(&disconnected) == BP_ERR_NO_PATH);

    BpDecomposition empty = {NULL, 0, NULL, 0};
    CHECK(bp_check_connected(&empty) == BP_ERR_ARGUMENT);
}

static double independently_calculated_length(const BpVisitPlan *plan,
        const BpDecomposition *decomposition, BpPoint start, BpPoint end) {
    BpPoint previous = start;
    double length = 0.0;
    for (size_t i = 0; i < plan->count; ++i) {
        const BpPath *path = &decomposition->cells[plan->order[i]].variants[plan->variants[i]];
        CHECK(path->count > 0);
        length += hypot(previous.x - path->points[0].x, previous.y - path->points[0].y);
        for (size_t point = 1; point < path->count; ++point) {
            length += hypot(path->points[point - 1].x - path->points[point].x,
                path->points[point - 1].y - path->points[point].y);
        }
        previous = path->points[path->count - 1];
    }
    length += hypot(previous.x - end.x, previous.y - end.y);
    return length;
}

static void check_plan_is_permutation(const BpVisitPlan *plan, size_t cell_count) {
    CHECK(plan->count == cell_count);
    int seen[4] = {0, 0, 0, 0};
    CHECK(cell_count <= sizeof(seen) / sizeof(seen[0]));
    for (size_t i = 0; i < plan->count; ++i) {
        CHECK(plan->order[i] < cell_count);
        CHECK(seen[plan->order[i]] == 0);
        seen[plan->order[i]] = 1;
        CHECK(plan->variants[i] < 4);
    }
    for (size_t i = 0; i < cell_count; ++i) CHECK(seen[i] == 1);
}

static void test_planner_cost_and_determinism(void) {
    BpDomain domain = rectangle_domain();
    BpCell cells[3] = {0};
    BpPoint cell_points[3][4] = {
        {{10.0, 10.0}, {10.0, 10.0}, {10.0, 10.0}, {10.0, 10.0}},
        {{40.0, 10.0}, {40.0, 10.0}, {40.0, 10.0}, {40.0, 10.0}},
        {{70.0, 10.0}, {70.0, 10.0}, {70.0, 10.0}, {70.0, 10.0}}
    };
    for (size_t cell = 0; cell < 3; ++cell) set_cell_paths(&cells[cell], cell_points[cell]);
    BpAdjacency edges[] = {{0, 1}, {1, 2}};
    BpDecomposition decomposition = {cells, 3, edges, 2};
    BpPoint start = {2.0, 4.0};
    BpPoint end = {95.0, 80.0};
    BpVisitPlan first = {0};
    BpVisitPlan second = {0};

    CHECK(bp_order_cells(&domain, &decomposition, start, end, &first) == BP_OK);
    CHECK(bp_order_cells(&domain, &decomposition, start, end, &second) == BP_OK);
    check_plan_is_permutation(&first, 3);
    CHECK(second.count == first.count);
    for (size_t i = 0; i < first.count; ++i) {
        CHECK(first.order[i] == second.order[i]);
        CHECK(first.variants[i] == second.variants[i]);
    }
    CHECK(first.length <= first.initial_length + 1e-9);
    NEAR(first.length, independently_calculated_length(&first, &decomposition, start, end), 1e-8);
    CHECK(first.length > 0.0);
    bp_free_plan(&first);
    bp_free_plan(&second);
}

static void test_each_variant_index(void) {
    BpDomain domain = rectangle_domain();
    const BpPoint start = {10.0, 50.0};
    const BpPoint end = {90.0, 50.0};

    for (size_t preferred = 0; preferred < 4; ++preferred) {
        BpCell cell = {0};
        BpPoint paths[4] = {
            {50.0, 100.0}, {50.0, 100.0}, {50.0, 100.0}, {50.0, 100.0}
        };
        paths[preferred] = (BpPoint){50.0, 50.0};
        set_cell_paths(&cell, paths);
        BpDecomposition decomposition = {&cell, 1, NULL, 0};
        BpVisitPlan plan = {0};

        CHECK(bp_order_cells(&domain, &decomposition, start, end, &plan) == BP_OK);
        CHECK(plan.count == 1);
        CHECK(plan.order[0] == 0);
        CHECK(plan.variants[0] == preferred);
        CHECK(plan.length <= plan.initial_length + 1e-9);
        NEAR(plan.length, independently_calculated_length(&plan, &decomposition, start, end), 1e-8);
        bp_free_plan(&plan);
    }
}

static void test_two_point_variant_cost(void) {
    BpDomain domain = rectangle_domain();
    BpCell cell = {0};
    BpPoint single_point_paths[3] = {
        {50.0, 100.0}, {50.0, 100.0}, {50.0, 100.0}
    };
    BpPoint preferred_path[] = {{30.0, 50.0}, {70.0, 50.0}};
    for (size_t variant = 0, single = 0; variant < 4; ++variant) {
        if (variant == 2) {
            cell.variants[variant].points = preferred_path;
            cell.variants[variant].count = 2;
        } else {
            cell.variants[variant].points = &single_point_paths[single++];
            cell.variants[variant].count = 1;
        }
    }
    BpDecomposition decomposition = {&cell, 1, NULL, 0};
    BpVisitPlan plan = {0};
    BpPoint start = {10.0, 50.0};
    BpPoint end = {90.0, 50.0};

    CHECK(bp_order_cells(&domain, &decomposition, start, end, &plan) == BP_OK);
    CHECK(plan.count == 1);
    CHECK(plan.variants[0] == 2);
    CHECK(plan.length <= plan.initial_length + 1e-9);
    NEAR(plan.length, independently_calculated_length(&plan, &decomposition, start, end), 1e-8);
    bp_free_plan(&plan);
}

static void test_empty_planner_input(void) {
    BpDomain domain = rectangle_domain();
    BpDecomposition decomposition = {NULL, 0, NULL, 0};
    BpVisitPlan plan = {0};
    CHECK(bp_order_cells(&domain, &decomposition,
        (BpPoint){1.0, 1.0}, (BpPoint){2.0, 2.0}, &plan) == BP_ERR_ARGUMENT);
}

int main(void) {
    test_connectivity();
    test_planner_cost_and_determinism();
    test_each_variant_index();
    test_two_point_variant_cost();
    test_empty_planner_input();
    return EXIT_SUCCESS;
}
