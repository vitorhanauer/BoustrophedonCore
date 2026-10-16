#include "test_support.h"

#include <float.h>

static const double EPS = 1e-8;

static int close_value(double a, double b) {
    return fabs(a - b) <= EPS;
}

static int same_point(BpPoint a, BpPoint b) {
    return close_value(a.x, b.x) && close_value(a.y, b.y);
}

static BpPoint point(double x, double y) {
    BpPoint result = {x, y};
    return result;
}

static double cross(BpPoint a, BpPoint b) {
    return a.x * b.y - a.y * b.x;
}

static double dot(BpPoint a, BpPoint b) {
    return a.x * b.x + a.y * b.y;
}

static int finite_path(const BpPath *path) {
    if (!path || !path->points || path->count < 2) return 0;
    for (size_t i = 0; i < path->count; ++i) {
        if (!isfinite(path->points[i].x) || !isfinite(path->points[i].y)) return 0;
    }
    return 1;
}

static void add_parameter(double *parameters, size_t *count, double value) {
    if (value < -EPS || value > 1.0 + EPS) return;
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;
    parameters[(*count)++] = value;
}

/* Split a candidate segment at every boundary crossing, then check each
 * resulting open interval. This catches excursions through a concavity while
 * allowing travel along the cell boundary. */
static int segment_in_cell(BpPolygon polygon, BpPoint a, BpPoint b) {
    double parameters[2 + 2 * 64];
    size_t parameter_count = 0;
    BpPoint direction = point(b.x - a.x, b.y - a.y);
    size_t vertex_count = bp_ring_count(polygon);
    if (vertex_count == 0 || vertex_count > 64) return 0;
    add_parameter(parameters, &parameter_count, 0.0);
    add_parameter(parameters, &parameter_count, 1.0);

    for (size_t i = 0; i < vertex_count; ++i) {
        BpPoint c = polygon.vertices[i];
        BpPoint d = polygon.vertices[(i + 1) % vertex_count];
        BpPoint edge = point(d.x - c.x, d.y - c.y);
        BpPoint offset = point(c.x - a.x, c.y - a.y);
        double denominator = cross(direction, edge);
        if (fabs(denominator) > EPS) {
            double t = cross(offset, edge) / denominator;
            double u = cross(offset, direction) / denominator;
            if (u >= -EPS && u <= 1.0 + EPS) add_parameter(parameters, &parameter_count, t);
        } else if (fabs(cross(offset, direction)) <= EPS) {
            double length_squared = dot(direction, direction);
            if (length_squared > EPS * EPS) {
                add_parameter(parameters, &parameter_count, dot(offset, direction) / length_squared);
                BpPoint d_offset = point(d.x - a.x, d.y - a.y);
                add_parameter(parameters, &parameter_count, dot(d_offset, direction) / length_squared);
            }
        }
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
    for (size_t i = 1; i < parameter_count; ++i) {
        if (parameters[i] - parameters[i - 1] <= EPS) continue;
        double t = (parameters[i] + parameters[i - 1]) * 0.5;
        BpPoint middle = point(a.x + t * direction.x, a.y + t * direction.y);
        if (bp_point_location(polygon, middle) < 0) return 0;
    }
    return bp_point_location(polygon, a) >= 0 && bp_point_location(polygon, b) >= 0;
}

static void check_path_inside(BpPolygon polygon, const BpPath *path) {
    CHECK(finite_path(path));
    for (size_t i = 1; i < path->count; ++i) {
        CHECK(segment_in_cell(polygon, path->points[i - 1], path->points[i]));
    }
}

/* A sweep advances across the cell once; a perimeter lap would backtrack. */
static void check_no_extra_contour(const BpPath *path, double span) {
    double horizontal_travel = 0.0;
    for (size_t i = 1; i < path->count; ++i)
        horizontal_travel += fabs(path->points[i].x - path->points[i - 1].x);
    CHECK(close_value(horizontal_travel, span));
}

static void check_reverse(const BpPath *path, const BpPath *reversed) {
    CHECK(path->count == reversed->count);
    for (size_t i = 0; i < path->count; ++i) {
        CHECK(same_point(path->points[i], reversed->points[path->count - i - 1]));
    }
}

static int point_on_path(BpPoint sample, const BpPath *path, double radius) {
    double radius_squared = radius * radius + EPS;
    for (size_t i = 1; i < path->count; ++i) {
        BpPoint start = path->points[i - 1];
        BpPoint vector = point(path->points[i].x - start.x, path->points[i].y - start.y);
        BpPoint offset = point(sample.x - start.x, sample.y - start.y);
        double length_squared = dot(vector, vector);
        double t = length_squared > EPS ? dot(offset, vector) / length_squared : 0.0;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        BpPoint nearest = point(start.x + t * vector.x, start.y + t * vector.y);
        double dx = sample.x - nearest.x;
        double dy = sample.y - nearest.y;
        if (dx * dx + dy * dy <= radius_squared) return 1;
    }
    return 0;
}

static void check_rectangle_coverage(BpPolygon polygon, const BpPath *path, double radius) {
    /* Analytic certificate: full-height vertical segments at x_i cover this
     * rectangle iff the extreme offsets <= radius and consecutive gaps <= 2r.
     * Then also check a grid, including the edges, as a separate regression. */
    double lanes[128]; size_t lane_count=0;
    for (size_t i=1;i<path->count;i++) {
        BpPoint a=path->points[i-1],b=path->points[i];
        if (a.x==b.x && fmin(a.y,b.y)<=EPS && fmax(a.y,b.y)>=4-EPS) {
            CHECK(lane_count<128);lanes[lane_count++]=a.x;
        }
    }
    CHECK(lane_count>0);
    for (size_t i=0;i<lane_count;i++) for (size_t j=i+1;j<lane_count;j++)
        if (lanes[j]<lanes[i]) { double t=lanes[i];lanes[i]=lanes[j];lanes[j]=t; }
    CHECK(lanes[0]<=radius+EPS);
    CHECK(6-lanes[lane_count-1]<=radius+EPS);
    for (size_t i=1;i<lane_count;i++) CHECK(lanes[i]-lanes[i-1]<=2*radius+EPS);
    for (int ix = 0; ix <= 60; ++ix) {
        for (int iy = 0; iy <= 40; ++iy) {
            BpPoint sample = point((double)ix / 10.0, (double)iy / 10.0);
            CHECK(bp_point_location(polygon, sample) >= 0);
            CHECK(point_on_path(sample, path, radius));
        }
    }
}

static void check_sweep_spacing(const BpPath *path, double xmin, double xmax, double step) {
    /* Every prescribed vertical cut must be represented in the path. */
    size_t cut_count = (size_t)ceil((xmax - xmin) / step) + 1;
    for (size_t cut = 0; cut < cut_count; ++cut) {
        double x = xmin + (double)cut * step;
        if (x > xmax) x = xmax;
        int found = 0;
        for (size_t i = 0; i < path->count; ++i) {
            if (close_value(path->points[i].x, x)) {
                found = 1;
                break;
            }
        }
        CHECK(found);
    }
    if (!close_value(xmin + (double)(cut_count - 1) * step, xmax)) {
        int found = 0;
        for (size_t i = 0; i < path->count; ++i) {
            if (close_value(path->points[i].x, xmax)) found = 1;
        }
        CHECK(found);
    }
}

static void exercise_cell(BpPolygon polygon, double step) {
    BpCell cell = {0};
    cell.boundary.points = (BpPoint *)malloc(polygon.count * sizeof(BpPoint));
    CHECK(cell.boundary.points != NULL);
    cell.boundary.count = polygon.count;
    for (size_t i = 0; i < polygon.count; ++i) cell.boundary.points[i] = polygon.vertices[i];

    CHECK(bp_cover_cell(&cell, step) == BP_OK);
    double xmin = polygon.vertices[0].x;
    double xmax = xmin;
    for (size_t i = 1; i < bp_ring_count(polygon); ++i) {
        if (polygon.vertices[i].x < xmin) xmin = polygon.vertices[i].x;
        if (polygon.vertices[i].x > xmax) xmax = polygon.vertices[i].x;
    }
    for (size_t variant = 0; variant < 4; ++variant) {
        check_path_inside(polygon, &cell.variants[variant]);
        check_no_extra_contour(&cell.variants[variant], xmax - xmin);
        check_sweep_spacing(&cell.variants[variant], xmin, xmax, step);
    }
    check_reverse(&cell.variants[0], &cell.variants[2]);
    check_reverse(&cell.variants[1], &cell.variants[3]);

    for (size_t variant = 0; variant < 4; ++variant) free(cell.variants[variant].points);
    free(cell.boundary.points);
}

static void test_rectangle(void) {
    const BpPoint vertices[] = {{0, 0}, {6, 0}, {6, 4}, {0, 4}};
    BpPolygon polygon = {vertices, sizeof(vertices) / sizeof(vertices[0])};
    BpCell cell = {0};
    cell.boundary.points = (BpPoint *)malloc(sizeof(vertices));
    CHECK(cell.boundary.points != NULL);
    for (size_t i = 0; i < polygon.count; ++i) cell.boundary.points[i] = vertices[i];
    cell.boundary.count = polygon.count;
    CHECK(bp_cover_cell(&cell, 1.0) == BP_OK);
    for (size_t variant = 0; variant < 4; ++variant) {
        check_path_inside(polygon, &cell.variants[variant]);
        check_no_extra_contour(&cell.variants[variant], 6.0);
        check_sweep_spacing(&cell.variants[variant], 0.0, 6.0, 1.0);
        check_rectangle_coverage(polygon, &cell.variants[variant], 0.5);
    }
    check_reverse(&cell.variants[0], &cell.variants[2]);
    check_reverse(&cell.variants[1], &cell.variants[3]);
    for (size_t variant = 0; variant < 4; ++variant) free(cell.variants[variant].points);
    free(cell.boundary.points);
}

static void test_invalid_steps(void) {
    const BpPoint vertices[] = {{0, 0}, {2, 0}, {2, 1}, {0, 1}};
    BpCell cell = {0};
    cell.boundary.points = (BpPoint *)malloc(sizeof(vertices));
    CHECK(cell.boundary.points != NULL);
    for (size_t i = 0; i < 4; ++i) cell.boundary.points[i] = vertices[i];
    cell.boundary.count = 4;
    CHECK(bp_cover_cell(&cell, 0.0) == BP_ERR_ARGUMENT);
    CHECK(bp_cover_cell(&cell, -1.0) == BP_ERR_ARGUMENT);
    CHECK(bp_cover_cell(&cell, NAN) == BP_ERR_ARGUMENT);
    CHECK(bp_cover_cell(&cell, INFINITY) == BP_ERR_ARGUMENT);
    free(cell.boundary.points);
}

int main(void) {
    const BpPoint triangle[] = {{0, 0}, {4, 0}, {2, 3}};
    const BpPoint concave[] = {{0, 0}, {5, 0}, {5, 4}, {3, 4}, {3, 2}, {2, 2}, {2, 4}, {0, 4}};
    const BpPoint narrow[] = {{1, 0}, {1.25, 0}, {1.25, 2}, {1, 2}};
    BpPolygon triangle_polygon = {triangle, sizeof(triangle) / sizeof(triangle[0])};
    BpPolygon concave_polygon = {concave, sizeof(concave) / sizeof(concave[0])};
    BpPolygon narrow_polygon = {narrow, sizeof(narrow) / sizeof(narrow[0])};

    test_rectangle();
    exercise_cell(triangle_polygon, 0.75);
    exercise_cell(concave_polygon, 0.8);
    exercise_cell(concave_polygon, 1.0);
    exercise_cell(narrow_polygon, 1.0);
    test_invalid_steps();
    return EXIT_SUCCESS;
}
