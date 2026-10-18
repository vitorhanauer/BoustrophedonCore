#include "test_support.h"
#include <stdint.h>

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

    BpPolygon oversized={square_vertices,SIZE_MAX};
    expect_api_error(BP_ERR_LIMIT,start,end,oversized,NULL,0,1,0);
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

/* Fractional intersections must survive decomposition, clipping and visibility,
 * including paths that follow an inclined boundary after rounding to double. */
static void test_inclined_pipeline(void) {
    const BpPoint triangle[]={{0,0},{7,1},{2,6}};
    const BpPoint general[]={ {5,0}, {-3.3632093666488725,5.825249499527395},
        {-1.6284242275864134,-2.8205134984557687} };
    const BpPoint concave[]={{0,0},{5,0},{5,4},{3,4},{3,2},{2,2},{2,4},{0,4}};
    BpPolygon fixtures[]={{triangle,3},{general,3},{concave,8}};
    for (size_t i=0;i<3;i++) {
        BpPolygon area=fixtures[i]; BpPath path={0};
        CHECK(generateWaypoints(area.vertices[0],area.vertices[1],area,NULL,0,0.7,25,&path)==BP_OK);
        check_safe_path(area,NULL,0,&path);
        freeWaypoints(&path);
    }
}

static void test_deterministic_general_polygons(void) {
    const double pi=acos(-1.0);
    for (size_t fixture=0;fixture<70;fixture++) {
        size_t n=3+fixture%7; double phase=(double)(fixture/7)*0.19;
        BpPoint vertices[9];
        for (size_t i=0;i<n;i++) {
            double angle=2*pi*(double)i/(double)n+phase;
            double radius=5+2*sin((double)i*2.1+phase);
            vertices[i]=(BpPoint){radius*cos(angle),radius*sin(angle)};
        }
        BpPolygon area={vertices,n}; BpPath path={0};
        CHECK(generateWaypoints(vertices[0],vertices[1],area,NULL,0,0.7,25,&path)==BP_OK);
        check_safe_path(area,NULL,0,&path);
        freeWaypoints(&path);
    }
}

static void test_general_obstacles_and_closed_route(void) {
    const BpPoint outer[]={{0,0},{20,0},{20,20},{0,20}};
    const BpPoint left[]={{4,8},{7,6},{9,9},{6,12}};
    const BpPoint right[]={{12,8},{15,7},{17,10},{14,13},{11,11}};
    BpPolygon area={outer,4},holes[]={{left,4},{right,5}};
    BpPoint start={1,1}; BpPath path={0};
    CHECK(generateWaypoints(start,start,area,holes,2,1.3,20,&path)==BP_OK);
    CHECK(bp_same_point(path.points[0],start));
    CHECK(bp_same_point(path.points[path.count-1],start));
    check_safe_path(area,holes,2,&path);freeWaypoints(&path);
    expect_api_error(BP_ERR_GEOMETRY,(BpPoint){6,9},start,area,holes,2,1,0);
    const BpPoint touching[]={{0,4},{3,4},{3,6},{0,6}};
    BpPolygon invalid={touching,4};
    expect_api_error(BP_ERR_GEOMETRY,start,start,area,&invalid,1,1,0);
}

/* The footprint must still reach edges and corners after moving the route. */
static double distance_to_route(BpPoint sample, const BpPath *path) {
    double best = INFINITY;
    for (size_t i = 1; i < path->count; ++i) {
        BpPoint a = path->points[i-1], b = path->points[i];
        double dx=b.x-a.x, dy=b.y-a.y, length=dx*dx+dy*dy;
        double t=length>0 ? ((sample.x-a.x)*dx+(sample.y-a.y)*dy)/length : 0;
        t=fmax(0,fmin(1,t));
        best=fmin(best,hypot(sample.x-a.x-t*dx,sample.y-a.y-t*dy));
    }
    return best;
}

static void test_public_boundary_clearance(void) {
    BpPolygon area={square_vertices,4}, obstacle={box_vertices,4};
    const double overlaps[]={0,25,75};
    for (size_t fixture=0;fixture<2;fixture++) {
        for (size_t k=0;k<3;k++) {
            BpPath path={0}; double width=1, overlap=overlaps[k];
            BpPoint start={1,1}, end={9,9};
            CHECK(generateWaypoints(start,end,area,fixture?&obstacle:NULL,
                fixture, width,overlap,&path)==BP_OK);
            check_safe_path(area,fixture?&obstacle:NULL,fixture,&path);
            for (size_t i=1;i+1<path.count;i++) {
                BpPoint p=path.points[i];
                CHECK(p.x>0 && p.x<10 && p.y>0 && p.y<10);
                if (fixture) CHECK(independent_location(obstacle,p)<0);
            }
            double margin=(1-overlap/100)/(2*sqrt(2.0));
            const BpPoint closed_box[]={{4,3},{6,3},{6,7},{4,7},{4,3}};
            BpPath box_boundary={(BpPoint *)closed_box,5};
            for (size_t i=0;i<path.count;i++) {
                BpPoint p=path.points[i];
                CHECK(fmin(fmin(p.x,10-p.x),fmin(p.y,10-p.y))>=margin-1e-8);
                if (fixture) CHECK(distance_to_route(p,&box_boundary)>=margin-1e-8);
            }
            if (fixture) {
                /* Checking both segment endpoint sets certifies separation
                 * of disjoint line segments, including diagonal connectors. */
                for (size_t i=0;i<4;i++)
                    CHECK(distance_to_route(closed_box[i],&path)>=margin-1e-8);
            }
            for (int x=0;x<=100;x++) for (int y=0;y<=100;y++) {
                BpPoint sample={(double)x/10,(double)y/10};
                if (fixture && independent_location(obstacle,sample)>0) continue;
                CHECK(distance_to_route(sample,&path)<=width/2+1e-8);
            }
            freeWaypoints(&path);
        }
    }
}

static void test_clearance_in_narrow_and_small_domains(void) {
    const BpPoint narrow[]={
        {0,0},{4,0},{4,1.9},{6,1.9},{6,0},{10,0},
        {10,4},{6,4},{6,2.1},{4,2.1},{4,4},{0,4}
    };
    const BpPoint tiny[]={{0,0},{0.1,0},{0.1,0.1},{0,0.1}};
    BpPolygon fixtures[]={{narrow,12},{tiny,4}};
    for (size_t i=0;i<2;i++) {
        BpPolygon area=fixtures[i];
        BpPoint start=i ? (BpPoint){0.02,0.02} : (BpPoint){1,2};
        BpPoint end=i ? (BpPoint){0.08,0.08} : (BpPoint){9,2};
        BpPath path={0};
        CHECK(generateWaypoints(start,end,area,NULL,0,1,0,&path)==BP_OK);
        check_safe_path(area,NULL,0,&path);
        CHECK(same_point(path.points[0],start));
        CHECK(same_point(path.points[path.count-1],end));
        for (size_t j=1;j+1<path.count;j++)
            CHECK(independent_location(area,path.points[j])==1);
        freeWaypoints(&path);
    }
}

static void test_boundary_terminals_use_short_approaches(void) {
    BpPolygon area={square_vertices,4}, obstacle={box_vertices,4};
    const BpPoint starts[]={{5,0},{4,5}};
    double margin=1/(2*sqrt(2.0));
    for (size_t fixture=0;fixture<2;fixture++) {
        BpPath path={0};
        CHECK(generateWaypoints(starts[fixture],(BpPoint){9,9},area,
            fixture ? &obstacle : NULL,fixture,1,0,&path)==BP_OK);
        CHECK(same_point(path.points[0],starts[fixture]));
        NEAR(hypot(path.points[1].x-starts[fixture].x,
            path.points[1].y-starts[fixture].y),margin,1e-8);
        check_safe_path(area,fixture ? &obstacle : NULL,fixture,&path);
        freeWaypoints(&path);
    }
}

static void test_overlap_spacing(void) {
    const BpPoint vertices[]={{0,0},{10,0},{10,10},{0,10}};
    double margin=1.5/(2*sqrt(2.0));
    double expected[]={margin,margin+1.5,margin+3,margin+4.5,
        margin+6,margin+7.5,10-margin};
    int seen[7]={0}; BpPath path={0};BpPolygon area={vertices,4};
    CHECK(generateWaypoints(vertices[0],vertices[2],area,NULL,0,2,25,&path)==BP_OK);
    for (size_t i=1;i<path.count;i++) {
        BpPoint a=path.points[i-1],b=path.points[i];
        if (a.x==b.x && fabs(fmin(a.y,b.y)-margin)<1e-12
                && fabs(fmax(a.y,b.y)-(10-margin))<1e-12) {
            size_t lane=0;
            while (lane<7 && fabs(a.x-expected[lane])>1e-12) lane++;
            CHECK(lane<7);seen[lane]=1;
        }
    }
    for (size_t i=0;i<7;i++) CHECK(seen[i]);
    freeWaypoints(&path);
}

static void test_sweep_finishes_without_extra_contour(void) {
    const BpPoint vertices[] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    BpPolygon area = {vertices, 4};
    BpPath route = {0};
    CHECK(generateWaypoints(vertices[0], vertices[2], area, NULL, 0,
        1.0, 0.0, &route) == BP_OK);
    check_safe_path(area, NULL, 0, &route);
    CHECK(same_point(route.points[0], vertices[0]));
    CHECK(same_point(route.points[route.count - 1], vertices[2]));
    /* Three inset cuts, two links and the two terminal approaches. */
    NEAR(independent_length(&route), 4*(2-1/sqrt(2.0))+1, 1e-9);
    freeWaypoints(&route);
}

static void test_closed_route_returns_after_last_cut(void) {
    BpPolygon area = {square_vertices, 4};
    BpPath route = {0};
    CHECK(generateWaypoints(square_vertices[0], square_vertices[0], area,
        NULL, 0, 1.0, 10.0, &route) == BP_OK);
    check_safe_path(area, NULL, 0, &route);
    CHECK(same_point(route.points[0], square_vertices[0]));
    CHECK(same_point(route.points[route.count - 1], square_vertices[0]));
    /* Bound sweep plus direct return; a perimeter lap exceeds this bound. */
    double margin=0.9/(2*sqrt(2.0)), span=10-2*margin;
    double gap=sqrt(1-4*margin*margin);
    double lanes=ceil(span/gap)+1;
    CHECK(independent_length(&route)<=(lanes+1)*span+sqrt(2.0)*span+0.9+1e-9);
    freeWaypoints(&route);
}

int main(void) {
    test_public_boundary_clearance();
    test_boundary_terminals_use_short_approaches();
    test_clearance_in_narrow_and_small_domains();
    test_closed_route_returns_after_last_cut();
    test_sweep_finishes_without_extra_contour();
    test_overlap_spacing();
    test_general_obstacles_and_closed_route();
    test_deterministic_general_polygons();
    test_inclined_pipeline();
    test_shortest_path_around_rectangle();
    test_shortest_path_through_concavity();
    test_stitch_connects_a_cell_path();
    test_public_pipeline_and_ring_normalization();
    test_public_api_rejects_invalid_inputs();
    return EXIT_SUCCESS;
}
