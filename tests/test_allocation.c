#include "test_support.h"
#include <stdint.h>

/* Linker fault injection leaves the production allocator and ABI unchanged. */
void *__real_malloc(size_t size);
void *__real_calloc(size_t count,size_t size);
void *__real_realloc(void *pointer,size_t size);
static size_t allocation_index,fail_at;
static int injecting;
static int fail_allocation(void) {
    return injecting && ++allocation_index==fail_at;
}
void *__wrap_malloc(size_t size) {
    return fail_allocation() ? NULL : __real_malloc(size);
}
void *__wrap_calloc(size_t count,size_t size) {
    return fail_allocation() ? NULL : __real_calloc(count,size);
}
void *__wrap_realloc(void *pointer,size_t size) {
    return fail_allocation() ? NULL : __real_realloc(pointer,size);
}
static BpStatus exercise(int connector,BpPath *output) {
    const BpPoint vertices[]={{0,0},{10,0},{10,10},{0,10}};
    const BpPoint block[]={{4,3},{6,3},{6,7},{4,7}};
    BpPolygon area={vertices,4},obstacle={block,4};
    BpDomain domain={area,&obstacle,1};
    if (connector==1) return bp_shortest_path(&domain,(BpPoint){2,5},(BpPoint){8,5},output,NULL);
    return generateWaypoints((BpPoint){0,0},(BpPoint){10,10},area,
        connector==2 ? &obstacle : NULL,connector==2 ? 1U : 0U,1,25,output);
}
int main(void) {
    for (int connector=0;connector<=2;connector++) {
        BpPath output={0};allocation_index=0;fail_at=SIZE_MAX;injecting=1;
        BpStatus status=exercise(connector,&output);
        injecting=0;CHECK(status==BP_OK);
        size_t allocation_count=allocation_index;CHECK(allocation_count>0);
        freeWaypoints(&output);
        for (size_t allocation=1;allocation<=allocation_count;allocation++) {
            allocation_index=0;fail_at=allocation;injecting=1;
            status=exercise(connector,&output);
            injecting=0;
            CHECK(status==BP_ERR_MEMORY);
            CHECK(output.points==NULL && output.count==0);
            freeWaypoints(&output);
        }
    }
    return EXIT_SUCCESS;
}
