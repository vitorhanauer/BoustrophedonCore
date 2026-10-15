#include "internal.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t *indices;
    size_t count;
} BoundaryChain;

static void free_chain(BoundaryChain *chain) {
    free(chain->indices);
    chain->indices = NULL;
    chain->count = 0;
}

/* Build one of the two boundary chains between the polygon's x extremes. */
static BpStatus make_chain(BpPolygon boundary, size_t first, size_t last,
    int direction, BoundaryChain *chain) {
    size_t n = bp_ring_count(boundary);
    if (n > SIZE_MAX / sizeof(*chain->indices)) return BP_ERR_LIMIT;
    chain->indices = malloc(n * sizeof(*chain->indices));
    if (!chain->indices) return BP_ERR_MEMORY;
    size_t index = first;
    for (;;) {
        if (chain->count >= n) {
            free_chain(chain);
            return BP_ERR_GEOMETRY;
        }
        chain->indices[chain->count++] = index;
        if (index == last) break;
        index = direction > 0 ? (index + 1) % n : (index + n - 1) % n;
    }
    return BP_OK;
}

/* Return the highest or lowest y where this chain meets a vertical cut. */
static int chain_y(const BoundaryChain *chain, BpPolygon boundary, double x,
    int highest, double *result) {
    int found = 0;
    double selected = 0.0;
    for (size_t i = 1; i < chain->count; ++i) {
        BpPoint a = boundary.vertices[chain->indices[i - 1]];
        BpPoint b = boundary.vertices[chain->indices[i]];
        double min_x = fmin(a.x, b.x), max_x = fmax(a.x, b.x);
        if (x < min_x || x > max_x) continue;
        double y;
        if (a.x == b.x) {
            y = highest ? fmax(a.y, b.y) : fmin(a.y, b.y);
        } else {
            long double fraction = ((long double)x - a.x) / ((long double)b.x - a.x);
            y = (double)((long double)a.y + fraction * ((long double)b.y - a.y));
        }
        if (!found || (highest ? y > selected : y < selected)) selected = y;
        found = 1;
    }
    if (found) *result = selected;
    return found;
}

/* Chain parameter, including vertical edges at a cut. Omitting vertices
 * whose x equals a cut would shortcut a vertical notch diagonally. */
static int chain_position(const BoundaryChain *chain, BpPolygon boundary,
    BpPoint p, long double *position) {
    for (size_t i=1;i<chain->count;i++) {
        BpPoint a=boundary.vertices[chain->indices[i-1]];
        BpPoint b=boundary.vertices[chain->indices[i]];
        if (p.x<fmin(a.x,b.x)||p.x>fmax(a.x,b.x)) continue;
        long double t;
        if (a.x==b.x) {
            if (p.y<fmin(a.y,b.y)||p.y>fmax(a.y,b.y)) continue;
            if (a.y==b.y) continue;
            t=((long double)p.y-a.y)/((long double)b.y-a.y);
        } else {
            t=((long double)p.x-a.x)/((long double)b.x-a.x);
            double y=(double)((long double)a.y+t*((long double)b.y-a.y));
            if (y!=p.y) continue;
        }
        *position=(long double)(i-1)+t;return 1;
    }
    return 0;
}
static BpStatus append_chain_interval(BpPath *path, const BoundaryChain *chain,
    BpPolygon boundary, double from_x, double to_x, int highest) {
    double from_y,to_y;
    if (!chain_y(chain,boundary,from_x,highest,&from_y)
        ||!chain_y(chain,boundary,to_x,highest,&to_y)) return BP_ERR_GEOMETRY;
    BpPoint first={from_x,from_y},last={to_x,to_y};
    long double first_pos,last_pos;
    if (!chain_position(chain,boundary,first,&first_pos)
        ||!chain_position(chain,boundary,last,&last_pos)||first_pos>last_pos) return BP_ERR_GEOMETRY;
    BpStatus status=bp_path_append(path,first);
    for (size_t i=0;i<chain->count && status==BP_OK;i++) {
        if ((long double)i>first_pos && (long double)i<last_pos)
            status=bp_path_append(path,boundary.vertices[chain->indices[i]]);
    }
    if (status==BP_OK) status=bp_path_append(path,last);
    return status;
}

static int slice_limits(BpPolygon boundary, double x, double *low, double *high) {
    size_t n = bp_ring_count(boundary);
    int found = 0;
    for (size_t i = 0; i < n; ++i) {
        BpPoint a = boundary.vertices[i], b = boundary.vertices[(i + 1) % n];
        if (x < fmin(a.x, b.x) || x > fmax(a.x, b.x)) continue;
        double y1, y2;
        if (a.x == b.x) {
            y1 = a.y;
            y2 = b.y;
        } else {
            long double fraction = ((long double)x - a.x) / ((long double)b.x - a.x);
            y1 = y2 = (double)((long double)a.y + fraction * ((long double)b.y - a.y));
        }
        double edge_low = fmin(y1, y2), edge_high = fmax(y1, y2);
        if (!found || edge_low < *low) *low = edge_low;
        if (!found || edge_high > *high) *high = edge_high;
        found = 1;
    }
    return found;
}

static int point_on_edge(BpPoint p, BpPoint a, BpPoint b) {
    long double dx = (long double)b.x - a.x, dy = (long double)b.y - a.y;
    long double px = (long double)p.x - a.x, py = (long double)p.y - a.y;
    long double scale=fabsl(dx)*(fabsl(p.y)+fabsl(a.y))
        +fabsl(dy)*(fabsl(p.x)+fabsl(a.x))
        +fabsl(px)*(fabsl(b.y)+fabsl(a.y))
        +fabsl(py)*(fabsl(b.x)+fabsl(a.x));
    long double tolerance=8*DBL_EPSILON*scale+16*LDBL_EPSILON*(fabsl(dx*py)+fabsl(dy*px));
    if (fabsl(dx*py-dy*px)>tolerance) return 0;
    long double ex=8*DBL_EPSILON*fmax(fabs(a.x),fabs(b.x));
    long double ey=8*DBL_EPSILON*fmax(fabs(a.y),fabs(b.y));
    return (long double)p.x>=fmin(a.x,b.x)-ex && (long double)p.x<=fmax(a.x,b.x)+ex
        && (long double)p.y>=fmin(a.y,b.y)-ey && (long double)p.y<=fmax(a.y,b.y)+ey;
}

/* Continue from the sweep endpoint around the complete boundary. */
static BpStatus append_boundary_contour(BpPath *path, BpPolygon boundary) {
    size_t n = bp_ring_count(boundary), edge = SIZE_MAX;
    BpPoint anchor = path->points[path->count - 1];
    for (size_t i = 0; i < n; ++i) {
        if (point_on_edge(anchor, boundary.vertices[i], boundary.vertices[(i + 1) % n])) {
            edge = i;
            break;
        }
    }
    if (edge == SIZE_MAX) return BP_ERR_GEOMETRY;
    BpStatus status;
    for (size_t offset = 1; offset <= n; ++offset) {
        status = bp_path_append(path, boundary.vertices[(edge + offset) % n]);
        if (status != BP_OK) return status;
    }
    return bp_path_append(path, anchor);
}

static void reverse_path(BpPath *path) {
    for (size_t i = 0, j = path->count ? path->count - 1 : 0; i < j; ++i, --j) {
        BpPoint temporary = path->points[i];
        path->points[i] = path->points[j];
        path->points[j] = temporary;
    }
}

BpStatus bp_cover_cell(BpCell *cell, double step) {
    if (!cell || !cell->boundary.points || cell->boundary.count < 3
        || !isfinite(step) || step <= 0.0) return BP_ERR_ARGUMENT;

    BpPolygon boundary = {cell->boundary.points, cell->boundary.count};
    size_t n = bp_ring_count(boundary);
    if (n < 3) return BP_ERR_GEOMETRY;
    double min_x = boundary.vertices[0].x, max_x = min_x;
    size_t min_index = 0, max_index = 0;
    for (size_t i = 0; i < n; ++i) {
        BpPoint p = boundary.vertices[i];
        if (!isfinite(p.x) || !isfinite(p.y)) return BP_ERR_GEOMETRY;
        if (p.x < min_x) { min_x = p.x; min_index = i; }
        if (p.x > max_x) { max_x = p.x; max_index = i; }
    }
    double span = max_x - min_x;
    if (!isfinite(span) || span <= 0.0) return BP_ERR_GEOMETRY;
    if (min_x + step <= min_x) return BP_ERR_LIMIT;
    long double ratio = (long double)span / (long double)step;
    if (!isfinite(ratio) || ratio > (long double)(SIZE_MAX - 1)) return BP_ERR_LIMIT;
    long double lane_count_ld = ceill(ratio) + 1.0L;
    if (lane_count_ld < 2.0L || lane_count_ld > (long double)(SIZE_MAX / sizeof(double))) return BP_ERR_LIMIT;
    size_t lane_count = (size_t)lane_count_ld;
    double *cuts = malloc(lane_count * sizeof(*cuts));
    if (!cuts) return BP_ERR_MEMORY;
    for (size_t i = 0; i < lane_count; ++i) {
        long double x = (long double)min_x + (long double)i * (long double)step;
        cuts[i] = i + 1 == lane_count || x >= (long double)max_x ? max_x : (double)x;
        if (!isfinite(cuts[i]) || (i && cuts[i] <= cuts[i - 1])) {
            free(cuts);
            return BP_ERR_LIMIT;
        }
    }

    BoundaryChain chains[2] = {{0}};
    BpStatus status = make_chain(boundary, min_index, max_index, 1, &chains[0]);
    if (status == BP_OK) status = make_chain(boundary, min_index, max_index, -1, &chains[1]);
    if (status != BP_OK) {
        free_chain(&chains[0]); free_chain(&chains[1]); free(cuts); return status;
    }
    double middle_x = (double)(((long double)min_x + max_x) / 2.0L);
    double chain_y0, chain_y1;
    if (!chain_y(&chains[0], boundary, middle_x, 1, &chain_y0)
        || !chain_y(&chains[1], boundary, middle_x, 1, &chain_y1)) {
        free_chain(&chains[0]); free_chain(&chains[1]); free(cuts); return BP_ERR_GEOMETRY;
    }
    int upper_chain = chain_y0 >= chain_y1 ? 0 : 1;
    int lower_chain = 1 - upper_chain;

    for (size_t variant = 0; variant < 2 && status == BP_OK; ++variant) {
        BpPath *path = &cell->variants[variant];
        int up_first = variant == 1;
        double low, high;
        if (!slice_limits(boundary, cuts[0], &low, &high)) { status = BP_ERR_GEOMETRY; break; }
        status = bp_path_append(path, up_first ? (BpPoint){cuts[0], high} : (BpPoint){cuts[0], low});
        if (status == BP_OK) status = bp_path_append(path, up_first ? (BpPoint){cuts[0], low} : (BpPoint){cuts[0], high});
        for (size_t i = 1; i < lane_count && status == BP_OK; ++i) {
            int connect_upper = ((i + (up_first ? 1U : 0U)) % 2U) == 1U;
            const BoundaryChain *chain = &chains[connect_upper ? upper_chain : lower_chain];
            status = append_chain_interval(path, chain, boundary, cuts[i - 1], cuts[i], connect_upper);
            if (status != BP_OK) break;
            if (!slice_limits(boundary, cuts[i], &low, &high)) { status = BP_ERR_GEOMETRY; break; }
            status = bp_path_append(path, connect_upper ? (BpPoint){cuts[i], low} : (BpPoint){cuts[i], high});
        }
        if (status == BP_OK) status = append_boundary_contour(path, boundary);
    }
    if (status == BP_OK) {
        cell->variants[2] = cell->variants[0];
        cell->variants[2].points = malloc(cell->variants[0].count * sizeof(BpPoint));
        cell->variants[3] = cell->variants[1];
        cell->variants[3].points = malloc(cell->variants[1].count * sizeof(BpPoint));
        if (!cell->variants[2].points || !cell->variants[3].points) status = BP_ERR_MEMORY;
        else {
            for (size_t i = 0; i < cell->variants[0].count; ++i) cell->variants[2].points[i] = cell->variants[0].points[i];
            for (size_t i = 0; i < cell->variants[1].count; ++i) cell->variants[3].points[i] = cell->variants[1].points[i];
            reverse_path(&cell->variants[2]);
            reverse_path(&cell->variants[3]);
        }
    }

    free_chain(&chains[0]); free_chain(&chains[1]); free(cuts);
    if (status != BP_OK) {
        for (size_t i = 0; i < 4; ++i) freeWaypoints(&cell->variants[i]);
    }
    return status;
}
