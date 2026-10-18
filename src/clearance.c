#include "internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* Unit edge normals point into free space, regardless of ring winding. */
static double free_side(BpPolygon ring, int obstacle) {
    size_t n=bp_ring_count(ring);
    long double area=0;
    BpPoint origin=ring.vertices[0];
    for (size_t i=1;i+1<n;i++)
        area+=bp_orientation(origin,ring.vertices[i],ring.vertices[i+1]);
    return (area>0 ? 1.0 : -1.0)*(obstacle ? -1.0 : 1.0);
}

static BpPoint offset_direction(BpPolygon ring, size_t index, double sign) {
    size_t n=bp_ring_count(ring);
    BpPoint a=ring.vertices[(index+n-1)%n], b=ring.vertices[index];
    BpPoint c=ring.vertices[(index+1)%n];
    double first=hypot(b.x-a.x,b.y-a.y), second=hypot(c.x-b.x,c.y-b.y);
    BpPoint u={sign*(a.y-b.y)/first,sign*(b.x-a.x)/first};
    BpPoint v={sign*(b.y-c.y)/second,sign*(c.x-b.x)/second};
    /* Solve dot(offset,u) = dot(offset,v) = 1. */
    double denominator=1+u.x*v.x+u.y*v.y;
    if (denominator<=0) return (BpPoint){INFINITY,INFINITY};
    return (BpPoint){(u.x+v.x)/denominator,(u.y+v.y)/denominator};
}

void bp_free_inset_domain(BpInsetDomain *domain) {
    freeWaypoints(&domain->area);
    if (domain->obstacles) {
        for (size_t i=0;i<domain->domain.obstacle_count;i++)
            free((void *)domain->obstacles[i].vertices);
    }
    free(domain->obstacles);
    *domain=(BpInsetDomain){0};
}

static BpStatus allocate_ring(BpPolygon source, BpPoint **vertices) {
    size_t n=bp_ring_count(source);
    if (n>SIZE_MAX/sizeof(**vertices)) return BP_ERR_LIMIT;
    *vertices=malloc(n*sizeof(**vertices));
    return *vertices ? BP_OK : BP_ERR_MEMORY;
}

/* Validate containment as well as simplicity: a valid offset ring alone can
 * jump across a narrow feature when the requested margin is too large. */
static BpStatus contained_ring(const BpDomain *source, BpPolygon ring, int *valid) {
    for (size_t i=0;i<ring.count;i++) {
        int visible=0;
        BpStatus status=bp_segment_visible(source,ring.vertices[i],
            ring.vertices[(i+1)%ring.count],&visible);
        if (status!=BP_OK) return status;
        if (!visible) { *valid=0; return BP_OK; }
    }
    return BP_OK;
}

BpStatus bp_inset_domain(const BpDomain *source, double spacing, BpInsetDomain *output) {
    double clearance=spacing/(2*sqrt(2.0));
    /* Limit displacement at sharp corners to half a lane spacing. */
    for (size_t r=0;r<=source->obstacle_count;r++) {
        BpPolygon ring=r ? source->obstacles[r-1] : source->area;
        double sign=free_side(ring,r!=0);
        for (size_t i=0;i<bp_ring_count(ring);i++) {
            BpPoint direction=offset_direction(ring,i,sign);
            double length=hypot(direction.x,direction.y);
            if (!isfinite(length)) return BP_ERR_LIMIT;
            clearance=fmin(clearance,(spacing/2)/length);
        }
    }
    BpInsetDomain inset={0};
    BpStatus status=allocate_ring(source->area,&inset.area.points);
    if (status!=BP_OK) return status;
    inset.area.count=bp_ring_count(source->area);
    inset.domain.area=(BpPolygon){inset.area.points,inset.area.count};
    inset.domain.obstacle_count=source->obstacle_count;
    if (source->obstacle_count) {
        if (source->obstacle_count>SIZE_MAX/sizeof(*inset.obstacles)) {
            status=BP_ERR_LIMIT; goto done;
        }
        inset.obstacles=calloc(source->obstacle_count,sizeof(*inset.obstacles));
        if (!inset.obstacles) { status=BP_ERR_MEMORY; goto done; }
        inset.domain.obstacles=inset.obstacles;
        for (size_t r=0;r<source->obstacle_count;r++) {
            BpPoint *points=NULL;
            status=allocate_ring(source->obstacles[r],&points);
            if (status!=BP_OK) goto done;
            inset.obstacles[r]=(BpPolygon){points,bp_ring_count(source->obstacles[r])};
        }
    }
    for (size_t attempt=0;attempt<64;attempt++,clearance/=2) {
        int valid=clearance>0;
        for (size_t r=0;r<=source->obstacle_count;r++) {
            BpPolygon ring=r ? source->obstacles[r-1] : source->area;
            BpPoint *points=r ? (BpPoint *)inset.obstacles[r-1].vertices : inset.area.points;
            double sign=free_side(ring,r!=0);
            for (size_t i=0;i<bp_ring_count(ring);i++) {
                BpPoint direction=offset_direction(ring,i,sign), p=ring.vertices[i];
                points[i]=(BpPoint){p.x+clearance*direction.x,p.y+clearance*direction.y};
                if (bp_same_point(points[i],p)) valid=0;
            }
        }
        if (!valid) { status=BP_ERR_LIMIT; goto done; }
        /* Crossing an opposite edge can invert a thin ring without making
         * it self-intersect. Reject offsets that reverse any source edge. */
        for (size_t r=0;r<=source->obstacle_count && valid;r++) {
            BpPolygon original=r ? source->obstacles[r-1] : source->area;
            BpPolygon moved=r ? inset.obstacles[r-1] : inset.domain.area;
            for (size_t i=0;i<moved.count;i++) {
                size_t j=(i+1)%moved.count;
                long double dx=(long double)original.vertices[j].x-original.vertices[i].x;
                long double dy=(long double)original.vertices[j].y-original.vertices[i].y;
                long double mx=(long double)moved.vertices[j].x-moved.vertices[i].x;
                long double my=(long double)moved.vertices[j].y-moved.vertices[i].y;
                if (dx*mx+dy*my<=0) valid=0;
            }
        }
        if (!valid) continue;
        status=bp_validate_domain(&inset.domain,inset.area.points[0],inset.area.points[0]);
        if (status==BP_ERR_MEMORY || status==BP_ERR_LIMIT) goto done;
        if (status!=BP_OK) continue;
        for (size_t r=0;r<=source->obstacle_count && valid;r++) {
            BpPolygon ring=r ? inset.obstacles[r-1] : inset.domain.area;
            status=contained_ring(source,ring,&valid);
            if (status!=BP_OK) goto done;
            if (r) {
                BpPolygon original=source->obstacles[r-1];
                for (size_t i=0;i<bp_ring_count(original);i++)
                    if (bp_point_location(ring,original.vertices[i])!=1) valid=0;
            }
        }
        if (!valid) continue;
        BpDecomposition check={0};
        status=bp_decompose(&inset.domain,&check);
        if (status==BP_OK) status=bp_check_connected(&check);
        bp_free_decomposition(&check);
        if (status==BP_ERR_MEMORY || status==BP_ERR_LIMIT) goto done;
        if (status!=BP_OK) continue;
        inset.clearance=clearance;
        *output=inset;
        return BP_OK;
    }
    status=BP_ERR_NO_PATH;
 done:
    bp_free_inset_domain(&inset);
    return status;
}

/* Explicit terminals may lie in the removed border strip. Project onto offset
 * edges, retaining visible vertices as alternatives if a projection is blocked. */
BpStatus bp_inset_terminal(const BpDomain *source, const BpDomain *inset,
        BpPoint terminal, BpPoint *anchor) {
    if (bp_point_free(inset,terminal)) { *anchor=terminal; return BP_OK; }
    double best=INFINITY;
    for (size_t r=0;r<=inset->obstacle_count;r++) {
        BpPolygon ring=r ? inset->obstacles[r-1] : inset->area;
        for (size_t i=0;i<bp_ring_count(ring);i++) {
            BpPoint a=ring.vertices[i], b=ring.vertices[(i+1)%bp_ring_count(ring)];
            long double dx=(long double)b.x-a.x, dy=(long double)b.y-a.y;
            long double t=(((long double)terminal.x-a.x)*dx
                +((long double)terminal.y-a.y)*dy)/(dx*dx+dy*dy);
            t=fmaxl(0,fminl(1,t));
            BpPoint candidates[]={a,{(double)((long double)a.x+t*dx),
                (double)((long double)a.y+t*dy)}};
            for (size_t candidate=0;candidate<2;candidate++) {
                BpPoint p=candidates[candidate];
                double distance=hypot(p.x-terminal.x,p.y-terminal.y);
                if (distance>=best || !bp_point_free(inset,p)) continue;
                int visible=0;
                BpStatus status=bp_segment_visible(source,terminal,p,&visible);
                if (status!=BP_OK) return status;
                if (visible) { best=distance; *anchor=p; }
            }
        }
    }
    return isfinite(best) ? BP_OK : BP_ERR_NO_PATH;
}
