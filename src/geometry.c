#include "internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

int bp_same_point(BpPoint a, BpPoint b) { return a.x == b.x && a.y == b.y; }
long double bp_orientation(BpPoint a, BpPoint b, BpPoint c) {
    return ((long double)b.x-a.x)*((long double)c.y-a.y)
         - ((long double)b.y-a.y)*((long double)c.x-a.x);
}
size_t bp_ring_count(BpPolygon ring) {
    if (ring.vertices && ring.count > 1 && bp_same_point(ring.vertices[0],ring.vertices[ring.count-1]))
        return ring.count-1;
    return ring.count;
}
static int orientation_sign(BpPoint a, BpPoint b, BpPoint c) {
    long double p=((long double)b.x-a.x)*((long double)c.y-a.y);
    long double q=((long double)b.y-a.y)*((long double)c.x-a.x);
    /* Account for rounded input/intersection coordinates, including points
     * near a vertex where product-only cancellation bounds vanish. */
    long double dx=(long double)b.x-a.x,dy=(long double)b.y-a.y;
    long double px=(long double)c.x-a.x,py=(long double)c.y-a.y;
    long double input_error=fabsl(dx)*(fabsl(c.y)+fabsl(a.y))
        +fabsl(dy)*(fabsl(c.x)+fabsl(a.x))
        +fabsl(px)*(fabsl(b.y)+fabsl(a.y))
        +fabsl(py)*(fabsl(b.x)+fabsl(a.x));
    long double value=p-q, tolerance=8*DBL_EPSILON*input_error
        +16*LDBL_EPSILON*(fabsl(p)+fabsl(q));
    return value > tolerance ? 1 : value < -tolerance ? -1 : 0;
}
static int on_segment(BpPoint a, BpPoint b, BpPoint p) {
    long double ex=8*DBL_EPSILON*fmax(fabs(a.x),fabs(b.x));
    long double ey=8*DBL_EPSILON*fmax(fabs(a.y),fabs(b.y));
    return orientation_sign(a,b,p)==0 && (long double)p.x>=fmin(a.x,b.x)-ex
        && (long double)p.x<=fmax(a.x,b.x)+ex
        && (long double)p.y>=fmin(a.y,b.y)-ey && (long double)p.y<=fmax(a.y,b.y)+ey;
}
static int segments_touch(BpPoint a, BpPoint b, BpPoint c, BpPoint d) {
    int ab_c=orientation_sign(a,b,c), ab_d=orientation_sign(a,b,d);
    int cd_a=orientation_sign(c,d,a), cd_b=orientation_sign(c,d,b);
    if (ab_c*ab_d<0 && cd_a*cd_b<0) return 1;
    return (!ab_c && on_segment(a,b,c)) || (!ab_d && on_segment(a,b,d))
        || (!cd_a && on_segment(c,d,a)) || (!cd_b && on_segment(c,d,b));
}
int bp_point_location(BpPolygon ring, BpPoint p) {
    size_t n=bp_ring_count(ring); int inside=0;
    for (size_t i=0;i<n;i++) {
        BpPoint a=ring.vertices[i], b=ring.vertices[(i+1)%n];
        if (on_segment(a,b,p)) return 0;
        if ((a.y>p.y)!=(b.y>p.y)) {
            int sign=orientation_sign(a,b,p);
            if ((b.y>a.y && sign>0) || (b.y<a.y && sign<0)) inside=!inside;
        }
    }
    return inside ? 1 : -1;
}
int bp_point_free(const BpDomain *domain, BpPoint p) {
    if (!isfinite(p.x) || !isfinite(p.y) || bp_point_location(domain->area,p)<0) return 0;
    for (size_t i=0;i<domain->obstacle_count;i++)
        if (bp_point_location(domain->obstacles[i],p)>0) return 0;
    return 1;
}
static BpStatus validate_ring(BpPolygon ring) {
    if (!ring.vertices || ring.count<3) return BP_ERR_ARGUMENT;
    if (ring.count>SIZE_MAX/sizeof(BpPoint)) return BP_ERR_LIMIT;
    size_t n=bp_ring_count(ring);
    if (n<3 || n>SIZE_MAX/sizeof(BpPoint)) return BP_ERR_GEOMETRY;
    long double area=0;
    BpPoint origin=ring.vertices[0];
    for (size_t i=0;i<n;i++) {
        BpPoint a=ring.vertices[i], b=ring.vertices[(i+1)%n];
        if (!isfinite(a.x)||!isfinite(a.y)) return BP_ERR_ARGUMENT;
        if (bp_same_point(a,b)) return BP_ERR_GEOMETRY;
        area+=bp_orientation(origin,a,b);
        BpPoint previous=ring.vertices[(i+n-1)%n];
        if (orientation_sign(previous,a,b)==0 && on_segment(previous,a,b)) return BP_ERR_GEOMETRY;
        for (size_t j=i+1;j<n;j++) {
            if (j==i+1 || (i==0 && j==n-1)) continue;
            if (segments_touch(a,b,ring.vertices[j],ring.vertices[(j+1)%n])) return BP_ERR_GEOMETRY;
        }
    }
    if (area==0 || !isfinite(area)) return BP_ERR_GEOMETRY;
    return BP_OK;
}
static int rings_touch(BpPolygon a, BpPolygon b) {
    size_t na=bp_ring_count(a), nb=bp_ring_count(b);
    for (size_t i=0;i<na;i++) for (size_t j=0;j<nb;j++)
        if (segments_touch(a.vertices[i],a.vertices[(i+1)%na],b.vertices[j],b.vertices[(j+1)%nb])) return 1;
    return 0;
}
BpStatus bp_validate_domain(const BpDomain *domain, BpPoint start, BpPoint end) {
    if (!domain || (domain->obstacle_count && !domain->obstacles)) return BP_ERR_ARGUMENT;
    if (domain->obstacle_count>SIZE_MAX/sizeof(BpPolygon)) return BP_ERR_LIMIT;
    BpStatus status=validate_ring(domain->area);
    if (status!=BP_OK) return status;
    for (size_t i=0;i<domain->obstacle_count;i++) {
        BpPolygon obstacle=domain->obstacles[i];
        status=validate_ring(obstacle);
        if (status!=BP_OK) return status;
        if (rings_touch(domain->area,obstacle) || bp_point_location(domain->area,obstacle.vertices[0])!=1) return BP_ERR_GEOMETRY;
        for (size_t j=0;j<i;j++) {
            BpPolygon other=domain->obstacles[j];
            if (rings_touch(obstacle,other) || bp_point_location(other,obstacle.vertices[0])>=0
                || bp_point_location(obstacle,other.vertices[0])>=0) return BP_ERR_GEOMETRY;
        }
    }
    if (!isfinite(start.x)||!isfinite(start.y)||!isfinite(end.x)||!isfinite(end.y)) return BP_ERR_ARGUMENT;
    if (!bp_point_free(domain,start)||!bp_point_free(domain,end)) return BP_ERR_GEOMETRY;
    double minx=domain->area.vertices[0].x,maxx=minx,miny=domain->area.vertices[0].y,maxy=miny;
    for (size_t i=1;i<bp_ring_count(domain->area);i++) {
        BpPoint p=domain->area.vertices[i]; minx=fmin(minx,p.x);maxx=fmax(maxx,p.x);miny=fmin(miny,p.y);maxy=fmax(maxy,p.y);
    }
    if (!isfinite(maxx-minx)||!isfinite(maxy-miny)||!isfinite(hypot(maxx-minx,maxy-miny))) return BP_ERR_LIMIT;
    return BP_OK;
}
static int compare_parameters(const void *left,const void *right) {
    long double a=*(const long double *)left,b=*(const long double *)right;
    return (a>b)-(a<b);
}
/* Split at every boundary intersection. Endpoints alone or a single midpoint
 * cannot establish visibility through a nonconvex polygon. */
BpStatus bp_segment_visible(const BpDomain *domain, BpPoint a, BpPoint b, int *visible) {
    if (!domain||!visible) return BP_ERR_ARGUMENT;
    *visible=0;
    if (!bp_point_free(domain,a)||!bp_point_free(domain,b)) return BP_OK;
    if (bp_same_point(a,b)) { *visible=1;return BP_OK; }
    size_t vertices=bp_ring_count(domain->area);
    for (size_t i=0;i<domain->obstacle_count;i++) {
        size_t n=bp_ring_count(domain->obstacles[i]);
        if (n>SIZE_MAX-vertices) return BP_ERR_LIMIT;
        vertices+=n;
    }
    if (vertices>(SIZE_MAX/sizeof(long double)-2)/2) return BP_ERR_LIMIT;
    long double *parameters=malloc((2*vertices+2)*sizeof(*parameters));
    if (!parameters) return BP_ERR_MEMORY;
    size_t count=0; parameters[count++]=0;parameters[count++]=1;
    long double rx=(long double)b.x-a.x,ry=(long double)b.y-a.y;
    for (size_t ring_index=0;ring_index<=domain->obstacle_count;ring_index++) {
        BpPolygon ring=ring_index ? domain->obstacles[ring_index-1] : domain->area;
        size_t n=bp_ring_count(ring);
        for (size_t i=0;i<n;i++) {
            BpPoint c=ring.vertices[i],d=ring.vertices[(i+1)%n];
            long double sx=(long double)d.x-c.x,sy=(long double)d.y-c.y;
            long double qx=(long double)c.x-a.x,qy=(long double)c.y-a.y;
            long double cross=rx*sy-ry*sx;
            if (cross!=0) {
                long double t=(qx*sy-qy*sx)/cross,u=(qx*ry-qy*rx)/cross;
                if (t>=0 && t<=1 && u>=0 && u<=1) parameters[count++]=t;
            } else if (qx*ry-qy*rx==0) {
                long double t1=fabsl(rx)>=fabsl(ry) ? qx/rx : qy/ry;
                long double t2=fabsl(rx)>=fabsl(ry) ? ((long double)d.x-a.x)/rx : ((long double)d.y-a.y)/ry;
                if (t1>=0&&t1<=1) parameters[count++]=t1;
                if (t2>=0&&t2<=1) parameters[count++]=t2;
            }
        }
    }
    qsort(parameters,count,sizeof(*parameters),compare_parameters);
    *visible=1;
    for (size_t i=1;i<count;i++) {
        if (parameters[i]==parameters[i-1]) continue;
        long double t=(parameters[i]+parameters[i-1])/2;
        BpPoint p={(double)((long double)a.x+t*rx),(double)((long double)a.y+t*ry)};
        if (!bp_point_free(domain,p)) { *visible=0;break; }
    }
    free(parameters);return BP_OK;
}
