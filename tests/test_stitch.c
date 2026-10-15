#include "test_support.h"

#include <float.h>

/* These helpers deliberately do not use bp_segment_visible or the production
 * point-location routine.  Each segment is split at every ring crossing and
 * every open interval is classified independently. */
static const BpPoint square_vertices[] = {
    {0, 0}, {10, 0}, {10, 10}, {0, 10}
};
static const BpPoint box_vertices[] = {
    {4, 3}, {6, 3}, {6, 7}, {4, 7}
};
static const BpPoint notch_vertices[] = {
    {0, 0}, {6, 0}, {6, 6}, {4, 6}, {4, 2}, {2, 2}, {2, 6}, {0, 6}
};

static BpPoint point(double x, double y) {
    return (BpPoint){x, y};
}

static size_t usable_count(BpPolygon ring) {
    size_t count = ring.count;
    if (count > 1 && fabs(ring.vertices[0].x - ring.vertices[count - 1].x) < 1e-10 &&
        fabs(ring.vertices[0].y - ring.vertices[count - 1].y) < 1e-10) {
        --count;
    }
    return count;
}

static double cross(BpPoint a, BpPoint b) {
    return a.x * b.y - a.y * b.x;
}

static double dot(BpPoint a, BpPoint b) {
    return a.x * b.x + a.y * b.y;
}

static int on_segment(BpPoint a, BpPoint b, BpPoint p) {
    const double scale = 1.0 + fabs(a.x) + fabs(a.y) + fabs(b.x) + fabs(b.y);
    BpPoint ab = point(b.x - a.x, b.y - a.y);
    BpPoint ap = point(p.x - a.x, p.y - a.y);
    if (fabs(cross(ab, ap)) > 1e-9 * scale) return 0;
    return dot(ap, ab) >= -1e-9 * scale && dot(ap, ab) <= dot(ab, ab) + 1e-9 * scale;
}

/* Returns -1 outside, 0 on the boundary, and 1 in the interior. */
static int independent_location(BpPolygon ring, BpPoint p) {
    int inside = 0;
    size_t count = usable_count(ring);
    for (size_t i = 0, j = count - 1; i < count; j = i++) {
        BpPoint a = ring.vertices[j];
        BpPoint b = ring.vertices[i];
        if (on_segment(a, b, p)) return 0;
        if ((a.y > p.y) != (b.y > p.y)) {
            double crossing_x = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if (crossing_x > p.x) inside = !inside;
        }
    }
    return inside ? 1 : -1;
}

static void add_parameter(double *values, size_t *count, double t) {
    if (t < -1e-9 || t > 1.0 + 1e-9) return;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    CHECK(*count < 1024);
    values[(*count)++] = t;
}

static void add_ring_intersections(BpPolygon ring, BpPoint a, BpPoint direction,
        double *parameters, size_t *parameter_count) {
    size_t count = usable_count(ring);
    for (size_t i = 0; i < count; ++i) {
        BpPoint c = ring.vertices[i];
        BpPoint d = ring.vertices[(i + 1) % count];
        BpPoint edge = point(d.x - c.x, d.y - c.y);
        BpPoint offset = point(c.x - a.x, c.y - a.y);
        double denominator = cross(direction, edge);
        if (fabs(denominator) > 1e-12) {
            double t = cross(offset, edge) / denominator;
            double u = cross(offset, direction) / denominator;
            if (u >= -1e-9 && u <= 1.0 + 1e-9) {
                add_parameter(parameters, parameter_count, t);
            }
        } else if (fabs(cross(offset, direction)) <= 1e-9) {
            double length_squared = dot(direction, direction);
            if (length_squared > 1e-24) {
                add_parameter(parameters, parameter_count,
                    dot(offset, direction) / length_squared);
                BpPoint d_offset = point(d.x - a.x, d.y - a.y);
                add_parameter(parameters, parameter_count,
                    dot(d_offset, direction) / length_squared);
            }
        }
    }
}

static int independently_safe_segment(BpPolygon area, const BpPolygon *obstacles,
        size_t obstacle_count, BpPoint a, BpPoint b) {
    double parameters[1024];
    size_t parameter_count = 0;
    BpPoint direction = point(b.x - a.x, b.y - a.y);
    add_parameter(parameters, &parameter_count, 0.0);
    add_parameter(parameters, &parameter_count, 1.0);
    add_ring_intersections(area, a, direction, parameters, &parameter_count);
    for (size_t i = 0; i < obstacle_count; ++i) {
        add_ring_intersections(obstacles[i], a, direction, parameters, &parameter_count);
    }
    for (size_t i = 0; i < parameter_count; ++i) {
        for (size_t j = i + 1; j < parameter_count; ++j) {
            if (parameters[j] < parameters[i]) {
                double temporary = parameters[i];
                parameters[i] = parameters[j];
                parameters[j] = temporary;
            }
        }
    }
    if (independent_location(area, a) < 0 || independent_location(area, b) < 0) return 0;
    for (size_t i = 1; i < parameter_count; ++i) {
        if (parameters[i] - parameters[i - 1] <= 1e-10) continue;
        double t = 0.5 * (parameters[i] + parameters[i - 1]);
        BpPoint middle = point(a.x + t * direction.x, a.y + t * direction.y);
        if (independent_location(area, middle) < 0) return 0;
        for (size_t obstacle = 0; obstacle < obstacle_count; ++obstacle) {
            /* Boundary travel is allowed; entering obstacle interior is not. */
            if (independent_location(obstacles[obstacle], middle) > 0) return 0;
        }
    }
    return 1;
}

static void check_safe_path(BpPolygon area, const BpPolygon *obstacles,
        size_t obstacle_count, const BpPath *path) {
    CHECK(path != NULL);
    CHECK(path->points != NULL);
    CHECK(path->count >= 2);
    for (size_t i = 0; i < path->count; ++i) {
        CHECK(isfinite(path->points[i].x) && isfinite(path->points[i].y));
        if (i > 0) {
            CHECK(independently_safe_segment(area, obstacles, obstacle_count,
                path->points[i - 1], path->points[i]));
        }
    }
}

static double independent_length(const BpPath *path) {
    double result = 0.0;
    for (size_t i = 1; i < path->count; ++i) {
        result += hypot(path->points[i].x - path->points[i - 1].x,
            path->points[i].y - path->points[i - 1].y);
    }
    return result;
}

static int same_point(BpPoint a, BpPoint b) {
    return fabs(a.x - b.x) <= 1e-9 && fabs(a.y - b.y) <= 1e-9;
}

static void test_shortest_path_around_rectangle(void) {
    BpPolygon area = {square_vertices, sizeof(square_vertices) / sizeof(square_vertices[0])};
    BpPolygon obstacle = {box_vertices, sizeof(box_vertices) / sizeof(box_vertices[0])};
    BpDomain domain = {area, &obstacle, 1};
    BpPoint start = point(2, 5);
    BpPoint end = point(8, 5);
    int visible = 1;
    BpPath route = {0};
    double distance = 0.0;

    CHECK(bp_segment_visible(&domain, start, end, &visible) == BP_OK);
    CHECK(visible == 0);
    CHECK(!independently_safe_segment(area, &obstacle, 1, start, end));
    CHECK(bp_shortest_path(&domain, start, end, &route, &distance) == BP_OK);
    CHECK(same_point(route.points[0], start));
    CHECK(same_point(route.points[route.count - 1], end));
    check_safe_path(area, &obstacle, 1, &route);
    NEAR(distance, 2.0 * sqrt(8.0) + 2.0, 1e-8);
    NEAR(independent_length(&route), distance, 1e-8);
    freeWaypoints(&route);

    /* A shortest route may touch and follow an obstacle boundary. */
    CHECK(bp_segment_visible(&domain, point(4, 3), point(6, 3), &visible) == BP_OK);
    CHECK(visible == 1);
    CHECK(independently_safe_segment(area, &obstacle, 1, point(4, 3), point(6, 3)));
}

static void test_shortest_path_through_concavity(void) {
    BpPolygon area = {notch_vertices, sizeof(notch_vertices) / sizeof(notch_vertices[0])};
    BpDomain domain = {area, NULL, 0};
    BpPoint start = point(1, 5);
    BpPoint end = point(5, 5);
    BpPath route = {0};
    double distance = 0.0;

    CHECK(!independently_safe_segment(area, NULL, 0, start, end));
    CHECK(bp_shortest_path(&domain, start, end, &route, &distance) == BP_OK);
    CHECK(same_point(route.points[0], start));
    CHECK(same_point(route.points[route.count - 1], end));
    check_safe_path(area, NULL, 0, &route);
    NEAR(distance, 2.0 * sqrt(10.0) + 2.0, 1e-8);
    NEAR(independent_length(&route), distance, 1e-8);
    freeWaypoints(&route);
}

static void test_stitch_connects_a_cell_path(void) {
    BpPolygon area = {square_vertices, sizeof(square_vertices) / sizeof(square_vertices[0])};
    BpDomain domain = {area, NULL, 0};
    BpPoint lane[] = {{3, 5}, {7, 5}};
    BpCell cell = {0};
    BpDecomposition decomposition = {&cell, 1, NULL, 0};
    size_t order[] = {0};
    size_t variants[] = {0};
    BpVisitPlan plan = {order, variants, 1, 0.0, 0.0};
    BpPath route = {0};
    BpPoint start = point(1, 1);
    BpPoint end = point(9, 9);

    cell.variants[0].points = lane;
    cell.variants[0].count = sizeof(lane) / sizeof(lane[0]);
    CHECK(bp_stitch(&domain, &decomposition, &plan, start, end, &route) == BP_OK);
    CHECK(same_point(route.points[0], start));
    CHECK(same_point(route.points[route.count - 1], end));
    check_safe_path(area, NULL, 0, &route);
    freeWaypoints(&route);
}

static void test_public_pipeline_and_ring_normalization(void) {
    const BpPoint closed_area[] = {
        {0, 0}, {10, 0}, {10, 10}, {0, 10}, {0, 0}
    };
    const BpPoint reversed_obstacle[] = {
        {4, 7}, {6, 7}, {6, 3}, {4, 3}, {4, 7}
    };
    BpPolygon area = {closed_area, sizeof(closed_area) / sizeof(closed_area[0])};
    BpPolygon obstacle = {reversed_obstacle,
        sizeof(reversed_obstacle) / sizeof(reversed_obstacle[0])};
    BpPoint start = point(1, 1);
    BpPoint end = point(9, 9);
    BpPath route = {0};

    CHECK(generateWaypoints(start, end, area, &obstacle, 1, 2.0, 0.0, &route) == BP_OK);
    CHECK(route.points != NULL && route.count > 2);
    CHECK(same_point(route.points[0], start));
    CHECK(same_point(route.points[route.count - 1], end));
    check_safe_path(area, &obstacle, 1, &route);
    freeWaypoints(&route);
    CHECK(route.points == NULL && route.count == 0);
    freeWaypoints(&route);
    CHECK(route.points == NULL && route.count == 0);
}

static void expect_api_error(BpStatus expected, BpPoint start, BpPoint end,
        BpPolygon area, const BpPolygon *obstacles, size_t obstacle_count,
        double width, double overlap) {
    BpPath output = {0};
    CHECK(generateWaypoints(start, end, area, obstacles, obstacle_count,
        width, overlap, &output) == expected);
    CHECK(output.points == NULL && output.count == 0);
}

static void test_public_api_rejects_invalid_inputs(void) {
    BpPolygon valid_area = {square_vertices,
        sizeof(square_vertices) / sizeof(square_vertices[0])};
    const BpPoint crossed_area_vertices[] = {
        {0, 0}, {10, 10}, {0, 10}, {10, 0}
    };
    const BpPoint invalid_obstacle_vertices[] = {
        {4, 3}, {6, 7}, {4, 7}, {6, 3}
    };
    BpPolygon crossed_area = {crossed_area_vertices,
        sizeof(crossed_area_vertices) / sizeof(crossed_area_vertices[0])};
    BpPolygon invalid_obstacle = {invalid_obstacle_vertices,
        sizeof(invalid_obstacle_vertices) / sizeof(invalid_obstacle_vertices[0])};
    BpPoint start = point(1, 1);
    BpPoint end = point(9, 9);

    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, 0.0, 0.0);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, -1.0, 0.0);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, NAN, 0.0);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, 1.0, 100.0);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, 1.0, -0.1);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 0, 1.0, NAN);
    expect_api_error(BP_ERR_ARGUMENT, start, end, valid_area, NULL, 1, 1.0, 0.0);
    CHECK(generateWaypoints(start, end, valid_area, NULL, 0, 1.0, 0.0, NULL) == BP_ERR_ARGUMENT);
    expect_api_error(BP_ERR_GEOMETRY, start, end, crossed_area, NULL, 0, 1.0, 0.0);
    expect_api_error(BP_ERR_GEOMETRY, start, end, valid_area, &invalid_obstacle, 1, 1.0, 0.0);
}

int main(void) {
    test_shortest_path_around_rectangle();
    test_shortest_path_through_concavity();
    test_stitch_connects_a_cell_path();
    test_public_pipeline_and_ring_normalization();
    test_public_api_rejects_invalid_inputs();
    return EXIT_SUCCESS;
}
