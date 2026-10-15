#ifndef BP_TEST_SUPPORT_H
#define BP_TEST_SUPPORT_H
#include "internal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); exit(EXIT_FAILURE); } } while (0)
#define NEAR(a,b,tol) CHECK(fabs((a)-(b)) <= (tol))
#endif
