#include <R.h>
#include <stdio.h>
#include <stdarg.h>
#include <functionWrappers.h>

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
};


void srand2(int seed) {
    return;
};
