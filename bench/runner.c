#define _POSIX_C_SOURCE 200809L

#include "boustrophedon.h"

#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

typedef enum { FIXTURE_SQUARE, FIXTURE_CONCAVE, FIXTURE_OBSTACLE } Fixture;

typedef struct {
    BpStatus status;
    size_t point_count;
    double path_distance;
    double elapsed_ms;
    long max_rss_kb;
} CaseResult;

static const char *const csv_header =
    "fixture,size,width,overlap,repeat,status,point_count,path_distance,elapsed_ms,max_rss_kb";

static void usage(FILE *stream, const char *program) {
    fprintf(stream,
        "Usage: %s [--all] [--size N --width W] [options]\n"
        "Run waypoint generation and report CSV to stdout.\n\n"
        "Options:\n"
        "  --all             Run the 14 standard size/width pairs (default mode)\n"
        "  --size N          Square side length (default: 100)\n"
        "  --width W         Sensor width (default: 1)\n"
        "  --overlap P       Overlap percentage in [0,100) (default: 0)\n"
        "  --fixture NAME    square, concave, or obstacle (default: square)\n"
        "  --repeat N        Repetitions per case (default: 1)\n"
        "  --no-header       Omit the CSV header\n"
        "  --help            Show this help\n",
        program);
}

static bool parse_double(const char *text, double *value) {
    char *end = NULL;
    errno = 0;
    double parsed = strtod(text, &end);
    if (errno || end == text || *end != '\0' || !isfinite(parsed)) return false;
    *value = parsed;
    return true;
}

static bool parse_size(const char *text, size_t *value) {
    if (!text[0] || text[0] == '-') return false;
    char *end = NULL;
    errno = 0;
    unsigned long long parsed = strtoull(text, &end, 10);
    if (errno || end == text || *end != '\0' || parsed == 0 || parsed > (unsigned long long)SIZE_MAX) return false;
    *value = (size_t)parsed;
    return true;
}

static bool parse_fixture(const char *text, Fixture *fixture) {
    if (strcmp(text, "square") == 0) *fixture = FIXTURE_SQUARE;
    else if (strcmp(text, "concave") == 0) *fixture = FIXTURE_CONCAVE;
    else if (strcmp(text, "obstacle") == 0) *fixture = FIXTURE_OBSTACLE;
    else return false;
    return true;
}

static const char *fixture_name(Fixture fixture) {
    switch (fixture) {
        case FIXTURE_SQUARE: return "square";
        case FIXTURE_CONCAVE: return "concave";
        case FIXTURE_OBSTACLE: return "obstacle";
    }
    return "square";
}

static bool write_all(int fd, const void *buffer, size_t size) {
    const unsigned char *bytes = buffer;
    size_t written = 0;
    while (written < size) {
        ssize_t result = write(fd, bytes + written, size - written);
        if (result < 0 && errno == EINTR) continue;
        if (result <= 0) return false;
        written += (size_t)result;
    }
    return true;
}

static bool read_all(int fd, void *buffer, size_t size) {
    unsigned char *bytes = buffer;
    size_t received = 0;
    while (received < size) {
        ssize_t result = read(fd, bytes + received, size - received);
        if (result < 0 && errno == EINTR) continue;
        if (result <= 0) return false;
        received += (size_t)result;
    }
    return true;
}

static CaseResult run_case(size_t size, double width, double overlap, Fixture fixture) {
    const double length = (double)size;
    BpPoint area_points[8];
    BpPolygon obstacles[1];
    BpPoint obstacle_points[4];
    size_t area_count = 4;
    area_points[0] = (BpPoint){0.0, 0.0};
    area_points[1] = (BpPoint){length, 0.0};
    area_points[2] = (BpPoint){length, length};
    area_points[3] = (BpPoint){0.0, length};

    if (fixture == FIXTURE_CONCAVE) {
        area_points[0] = (BpPoint){0.0, 0.0};
        area_points[1] = (BpPoint){length, 0.0};
        area_points[2] = (BpPoint){length, length};
        area_points[3] = (BpPoint){0.65 * length, length};
        area_points[4] = (BpPoint){0.65 * length, 0.65 * length};
        area_points[5] = (BpPoint){0.35 * length, 0.65 * length};
        area_points[6] = (BpPoint){0.35 * length, length};
        area_points[7] = (BpPoint){0.0, length};
        area_count = 8;
    } else if (fixture == FIXTURE_OBSTACLE) {
        obstacle_points[0] = (BpPoint){0.4 * length, 0.4 * length};
        obstacle_points[1] = (BpPoint){0.6 * length, 0.4 * length};
        obstacle_points[2] = (BpPoint){0.6 * length, 0.6 * length};
        obstacle_points[3] = (BpPoint){0.4 * length, 0.6 * length};
        obstacles[0] = (BpPolygon){obstacle_points, 4};
    }

    const BpPolygon area = {area_points, area_count};
    const BpPoint start = {0.0, 0.0};
    const BpPoint end = {length, length};
    const BpPolygon *obstacle_arg = fixture == FIXTURE_OBSTACLE ? obstacles : NULL;
    const size_t obstacle_count = fixture == FIXTURE_OBSTACLE ? 1U : 0U;
    BpPath path = {0};
    struct timespec before, after;
    (void)clock_gettime(CLOCK_MONOTONIC, &before);
    CaseResult result = {0};
    result.status = generateWaypoints(start, end, area, obstacle_arg,
                                      obstacle_count, width, overlap, &path);
    (void)clock_gettime(CLOCK_MONOTONIC, &after);
    result.elapsed_ms = ((double)(after.tv_sec - before.tv_sec) * 1000.0) +
                        ((double)(after.tv_nsec - before.tv_nsec) / 1000000.0);

    if (result.status == BP_OK) {
        result.point_count = path.count;
        bool valid = path.count >= 2;
        for (size_t i = 0; valid && i < path.count; ++i) {
            const BpPoint point = path.points[i];
            valid = isfinite(point.x) && isfinite(point.y) &&
                    point.x >= 0.0 && point.x <= length &&
                    point.y >= 0.0 && point.y <= length;
        }
        if (valid) {
            const double epsilon = fmax(1.0, length) * 1e-9;
            valid = hypot(path.points[0].x - start.x, path.points[0].y - start.y) <= epsilon &&
                    hypot(path.points[path.count - 1].x - end.x,
                          path.points[path.count - 1].y - end.y) <= epsilon;
        }
        if (!valid) result.status = BP_ERR_GEOMETRY;
        else {
            for (size_t i = 1; i < path.count; ++i)
                result.path_distance += hypot(path.points[i].x - path.points[i - 1].x,
                                              path.points[i].y - path.points[i - 1].y);
            if (!isfinite(result.path_distance)) result.status = BP_ERR_GEOMETRY;
        }
    }
    freeWaypoints(&path);

    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) result.max_rss_kb = usage.ru_maxrss;
    return result;
}

static const char *status_name(BpStatus status) {
    switch (status) {
        case BP_OK: return "success";
        case BP_ERR_ARGUMENT: return "invalid_argument";
        case BP_ERR_GEOMETRY: return "invalid_geometry";
        case BP_ERR_NO_PATH: return "no_path";
        case BP_ERR_LIMIT: return "limit_exceeded";
        case BP_ERR_MEMORY: return "allocation_failed";
    }
    return "unknown_status";
}

static void print_row(Fixture fixture, size_t size, double width, double overlap,
                      size_t repetition, const CaseResult *result) {
    printf("%s,%zu,%.17g,%.17g,%zu,%s,%zu,%.17g,%.6f,%ld\n",
           fixture_name(fixture), size, width, overlap, repetition,
           status_name(result->status), result->point_count,
           result->path_distance, result->elapsed_ms, result->max_rss_kb);
}

static bool run_isolated(Fixture fixture, size_t size, double width, double overlap,
                         CaseResult *result) {
    int pipe_fds[2];
    if (pipe(pipe_fds) != 0) return false;
    pid_t child = fork();
    if (child < 0) {
        close(pipe_fds[0]); close(pipe_fds[1]);
        return false;
    }
    if (child == 0) {
        close(pipe_fds[0]);
        CaseResult child_result = run_case(size, width, overlap, fixture);
        bool wrote = write_all(pipe_fds[1], &child_result, sizeof(child_result));
        close(pipe_fds[1]);
        _exit(wrote ? 0 : 2);
    }
    close(pipe_fds[1]);
    bool received = read_all(pipe_fds[0], result, sizeof(*result));
    close(pipe_fds[0]);
    int wait_status = 0;
    while (waitpid(child, &wait_status, 0) < 0) {
        if (errno != EINTR) return false;
    }
    return received && WIFEXITED(wait_status) && WEXITSTATUS(wait_status) == 0;
}

int main(int argc, char **argv) {
    bool all = false, no_header = false;
    bool size_set = false, width_set = false;
    size_t size = 100, repetitions = 1;
    double width = 1.0, overlap = 0.0;
    Fixture fixture = FIXTURE_SQUARE;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) { usage(stdout, argv[0]); return 0; }
        if (strcmp(argv[i], "--no-header") == 0) { no_header = true; continue; }
        if (strcmp(argv[i], "--all") == 0) { all = true; continue; }
        if (i + 1 >= argc) { usage(stderr, argv[0]); return 2; }
        const char *value = argv[++i];
        if (strcmp(argv[i - 1], "--size") == 0) {
            if (!parse_size(value, &size)) { fprintf(stderr, "Invalid --size: %s\n", value); return 2; }
            size_set = true;
        } else if (strcmp(argv[i - 1], "--width") == 0) {
            if (!parse_double(value, &width) || width <= 0.0) { fprintf(stderr, "Invalid --width: %s\n", value); return 2; }
            width_set = true;
        } else if (strcmp(argv[i - 1], "--overlap") == 0) {
            if (!parse_double(value, &overlap) || overlap < 0.0 || overlap >= 100.0) { fprintf(stderr, "Invalid --overlap: %s\n", value); return 2; }
        } else if (strcmp(argv[i - 1], "--repeat") == 0) {
            if (!parse_size(value, &repetitions)) { fprintf(stderr, "Invalid --repeat: %s\n", value); return 2; }
        } else if (strcmp(argv[i - 1], "--fixture") == 0) {
            if (!parse_fixture(value, &fixture)) { fprintf(stderr, "Invalid --fixture: %s\n", value); return 2; }
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i - 1]);
            usage(stderr, argv[0]);
            return 2;
        }
    }
    if (all && (size_set || width_set)) {
        fprintf(stderr, "--all cannot be combined with --size or --width\n");
        return 2;
    }

    if (!size_set && !width_set) all=true;
    if (!no_header) puts(csv_header);
    bool failed = false;
    if (all) {
        const size_t sizes[] = {10, 100, 1000, 10000};
        const double widths[] = {1.0, 5.0, 10.0, 50.0};
        for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
            const size_t width_count = si == 0 ? 2U : 4U;
            for (size_t wi = 0; wi < width_count; ++wi) {
                const double case_width = widths[wi];
                for (size_t repetition = 0; repetition < repetitions; ++repetition) {
                    CaseResult result = {BP_ERR_MEMORY, 0, 0.0, 0.0, 0};
                    if (!run_isolated(fixture, sizes[si], case_width, overlap, &result)) {
                        result.status = BP_ERR_MEMORY;
                        fprintf(stderr, "Child failed for %s size=%zu width=%.17g repeat=%zu\n",
                                fixture_name(fixture), sizes[si], case_width, repetition+1);
                    }
                    print_row(fixture, sizes[si], case_width, overlap, repetition+1, &result);
                    if (result.status != BP_OK) failed = true;
                }
            }
        }
    } else {
        for (size_t repetition = 0; repetition < repetitions; ++repetition) {
            CaseResult result = {BP_ERR_MEMORY, 0, 0.0, 0.0, 0};
            if (!run_isolated(fixture,size,width,overlap,&result)) {
                result.status=BP_ERR_MEMORY;
                fprintf(stderr,"Child failed for custom case\n");
            }
            print_row(fixture, size, width, overlap, repetition+1, &result);
            if (result.status != BP_OK) failed = true;
        }
    }
    return failed ? 1 : 0;
}
