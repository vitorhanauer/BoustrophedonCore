#include "boustrophedon.h"

#include <math.h>
#include <stdio.h>

int main(void) {
    const BpPoint vertices[] = {
        {0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}
    };
    const BpPolygon area = {vertices, sizeof(vertices) / sizeof(vertices[0])};
    const BpPoint start = {0.0, 0.0};
    const BpPoint end = {10.0, 10.0};
    BpPath path = {0};
    const BpStatus status = generateWaypoints(
        start, end, area, NULL, 0, 1.5, 0.0, &path);

    if (status != BP_OK) {
        fprintf(stderr, "Failed to generate waypoints: %s\n", bpStatusString(status));
        return 1;
    }

    double distance = 0.0;
    for (size_t i = 1; i < path.count; ++i) {
        distance += hypot(path.points[i].x - path.points[i - 1].x,
                          path.points[i].y - path.points[i - 1].y);
    }
    printf("Generated %zu waypoints; path distance %.6f\n", path.count, distance);
    freeWaypoints(&path);
    return 0;
}
