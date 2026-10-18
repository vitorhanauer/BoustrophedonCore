#ifndef BP_INTERNAL_H
#define BP_INTERNAL_H
#include "boustrophedon.h"
#include <stddef.h>

typedef struct {
    BpPolygon area;
    const BpPolygon *obstacles;
    size_t obstacle_count;
} BpDomain;
typedef struct {
    BpDomain domain;
    BpPath area;
    BpPolygon *obstacles;
    double clearance;
} BpInsetDomain;
BpStatus bp_inset_domain(const BpDomain *source, double spacing, BpInsetDomain *output);
void bp_free_inset_domain(BpInsetDomain *domain);
BpStatus bp_inset_terminal(const BpDomain *source, const BpDomain *inset,
    BpPoint terminal, BpPoint *anchor);
typedef struct {
    BpPath boundary; /* CCW or CW simple x-monotone polygon, owned. */
    BpPath variants[4]; /* two sweep senses, then their reversals. */
} BpCell;
typedef struct { size_t from, to; } BpAdjacency;
typedef struct {
    BpCell *cells;
    size_t count;
    BpAdjacency *adjacencies;
    size_t adjacency_count;
} BpDecomposition;
typedef struct {
    size_t *order;
    size_t *variants;
    size_t count;
    double initial_length;
    double length;
} BpVisitPlan;

int bp_same_point(BpPoint a, BpPoint b);
long double bp_orientation(BpPoint a, BpPoint b, BpPoint c);
size_t bp_ring_count(BpPolygon ring);
int bp_point_location(BpPolygon ring, BpPoint point); /* -1 exterior, 0 boundary, 1 interior */
int bp_point_free(const BpDomain *domain, BpPoint point);
BpStatus bp_validate_domain(const BpDomain *domain, BpPoint start, BpPoint end);
BpStatus bp_segment_visible(const BpDomain *domain, BpPoint a, BpPoint b, int *visible);
BpStatus bp_path_append(BpPath *path, BpPoint point);
BpStatus bp_path_extend(BpPath *path, const BpPath *other);
double bp_path_length(const BpPath *path);
void bp_free_decomposition(BpDecomposition *decomposition);
void bp_free_plan(BpVisitPlan *plan);

BpStatus bp_decompose(const BpDomain *domain, BpDecomposition *output);
BpStatus bp_cover_cell(BpCell *cell, double step);
BpStatus bp_check_connected(const BpDecomposition *decomposition);
BpStatus bp_order_cells(const BpDomain *domain, const BpDecomposition *decomposition,
    BpPoint start, BpPoint end, BpVisitPlan *output);
/* Shortest safe connector, allocated output; NULL output requests only distance. */
BpStatus bp_shortest_path(const BpDomain *domain, BpPoint start, BpPoint end,
    BpPath *output, double *distance);
BpStatus bp_stitch(const BpDomain *domain, const BpDecomposition *decomposition,
    const BpVisitPlan *plan, BpPoint start, BpPoint end, BpPath *output);
#endif
