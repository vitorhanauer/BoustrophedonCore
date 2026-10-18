#include "internal.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

void freeWaypoints(BpPath *path) {
    if (path) { free(path->points);path->points=NULL;path->count=0; }
}
const char *bpStatusString(BpStatus status) {
    switch(status) {
        case BP_OK:return "success";
        case BP_ERR_ARGUMENT:return "invalid argument";
        case BP_ERR_GEOMETRY:return "invalid geometry";
        case BP_ERR_NO_PATH:return "no path";
        case BP_ERR_LIMIT:return "representation limit exceeded";
        case BP_ERR_MEMORY:return "allocation failed";
    }
    return "unknown status";
}
BpStatus bp_path_append(BpPath *path, BpPoint point) {
    if (!path || !isfinite(point.x)||!isfinite(point.y)) return BP_ERR_ARGUMENT;
    if (path->count && bp_same_point(path->points[path->count-1],point)) return BP_OK;
    if (path->count>=SIZE_MAX/sizeof(*path->points)) return BP_ERR_LIMIT;
    BpPoint *points=realloc(path->points,(path->count+1)*sizeof(*points));
    if (!points) return BP_ERR_MEMORY;
    path->points=points;path->points[path->count++]=point;return BP_OK;
}
BpStatus bp_path_extend(BpPath *path, const BpPath *other) {
    for (size_t i=0;i<other->count;i++) {
        BpStatus status=bp_path_append(path,other->points[i]);
        if (status!=BP_OK) return status;
    }
    return BP_OK;
}
double bp_path_length(const BpPath *path) {
    double length=0;
    for (size_t i=1;i<path->count;i++) length+=hypot(path->points[i].x-path->points[i-1].x,path->points[i].y-path->points[i-1].y);
    return length;
}
void bp_free_decomposition(BpDecomposition *decomposition) {
    if (!decomposition) return;
    for (size_t i=0;i<decomposition->count;i++) {
        freeWaypoints(&decomposition->cells[i].boundary);
        for (size_t j=0;j<4;j++) freeWaypoints(&decomposition->cells[i].variants[j]);
    }
    free(decomposition->cells);free(decomposition->adjacencies);
    *decomposition=(BpDecomposition){0};
}
void bp_free_plan(BpVisitPlan *plan) {
    if (!plan) return;
    free(plan->order);free(plan->variants);*plan=(BpVisitPlan){0};
}
BpStatus generateWaypoints(BpPoint start, BpPoint end, BpPolygon area,
    const BpPolygon *obstacles,size_t obstacle_count,double width,double overlap_percent,BpPath *output) {
    if (!output) return BP_ERR_ARGUMENT;
    *output=(BpPath){0};
    if (!isfinite(width)||width<=0||!isfinite(overlap_percent)||overlap_percent<0||overlap_percent>=100) return BP_ERR_ARGUMENT;
    double step=width*(1-overlap_percent/100);
    if (!isfinite(step)||step<=0) return BP_ERR_LIMIT;
    BpDomain domain={area,obstacles,obstacle_count};
    BpStatus status=bp_validate_domain(&domain,start,end);
    if (status!=BP_OK) return status;
    BpDecomposition decomposition={0};BpVisitPlan plan={0};BpPath result={0};
    BpInsetDomain inset={0};
    BpPoint inset_start, inset_end;
    status=bp_inset_domain(&domain,step,&inset);
    if (status!=BP_OK) goto done;
    status=bp_inset_terminal(&domain,&inset.domain,start,&inset_start);
    if (status!=BP_OK) goto done;
    status=bp_inset_terminal(&domain,&inset.domain,end,&inset_end);
    if (status!=BP_OK) goto done;
    /* At a horizontal edge, the most distant sample lies midway between
     * lanes. Keep its distance sqrt((gap/2)^2 + clearance^2) <= width/2. */
    double relative_clearance=2*inset.clearance/width;
    step=fmin(step,width*sqrt(1-relative_clearance*relative_clearance));
    status=bp_decompose(&inset.domain,&decomposition);
    if (status!=BP_OK) goto done;
    for (size_t i=0;i<decomposition.count;i++) {
        status=bp_cover_cell(&decomposition.cells[i],step);
        if (status!=BP_OK) goto done;
    }
    status=bp_order_cells(&inset.domain,&decomposition,inset_start,inset_end,&plan);
    if (status!=BP_OK) goto done;
    BpPath coverage={0};
    status=bp_stitch(&inset.domain,&decomposition,&plan,inset_start,inset_end,&coverage);
    if (status==BP_OK) status=bp_path_append(&result,start);
    if (status==BP_OK) status=bp_path_extend(&result,&coverage);
    if (status==BP_OK) status=bp_path_append(&result,end);
    freeWaypoints(&coverage);
    if (status!=BP_OK) goto done;
    if (!isfinite(bp_path_length(&result))) { status=BP_ERR_LIMIT;goto done; }
    *output=result;result=(BpPath){0};
 done:
    freeWaypoints(&result);bp_free_plan(&plan);bp_free_decomposition(&decomposition);
    bp_free_inset_domain(&inset);
    return status;
}
