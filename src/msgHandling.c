#include <R.h>
#include <stdio.h>
#include <stdarg.h>

int fprintf2(FILE *stream, const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    if (stream == stdout)
        Rprintf(format, argptr);
    else if (stream == stderr) 
        REprintf(format, argptr);
    else
        vfprintf(stream, format, argptr);
    va_end(argptr);
}

int printf2(const char *format, ...) {
    va_list argptr;
    va_start(argptr, format);
    Rprintf(format, argptr);
    va_end(argptr);
    
}