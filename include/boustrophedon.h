#ifndef BOUSTROPHEDON_H
#define BOUSTROPHEDON_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double x, y; } BpPoint;
typedef struct { const BpPoint *vertices; size_t count; } BpPolygon;
typedef struct { BpPoint *points; size_t count; } BpPath;
typedef enum {
    BP_OK = 0, BP_ERR_ARGUMENT, BP_ERR_GEOMETRY, BP_ERR_NO_PATH,
    BP_ERR_LIMIT, BP_ERR_MEMORY
} BpStatus;

/* Inputs are borrowed. output must be empty; success transfers ownership.
 * overlap_percent is in [0,100). Coverage and connectors are offset from
 * original boundaries using the lane spacing; see README for the margin model.
 * Explicit start/end points are preserved and may lie on original boundaries.
 * Release successful results with freeWaypoints before reusing output. */
BpStatus generateWaypoints(BpPoint start, BpPoint end, BpPolygon area,
    const BpPolygon *obstacles, size_t obstacle_count, double width,
    double overlap_percent, BpPath *output);
void freeWaypoints(BpPath *path);
const char *bpStatusString(BpStatus status);
#ifdef __cplusplus
}
#endif
#endif
