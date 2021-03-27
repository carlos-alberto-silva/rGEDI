#include <R.h>
#include <stdio.h>
#include <stdarg.h>
#include "tools/functionWrappers.h"
#include <inttypes.h>

int msgf(const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    Rvprintf(format, argptr);
    va_end(argptr);
    return(0);
}

int errorf(const char *format, ...) {
    va_list argptr;
    va_start(argptr, format);
    REvprintf(format, argptr);
    va_end(argptr);
    return(0);
}

float frand() {
    return((float)unif_rand());
}

int rand2() {
    return((int)(unif_rand() * INT32_MAX));
}

void srand2(int seed) {
    return;
}
