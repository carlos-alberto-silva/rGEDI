#include <stdlib.h>
#include <string.h>
#include <Rinternals.h>
#include "argParse.h"

void allocAndCopy(char **argv, int *pos, const char *value)
{
    size_t len = strlen(value) + 1;
    argv[*pos] = malloc(len * sizeof(*value));
    strcpy(argv[(*pos)++], value);
}

void putArgName(char **argv, int *pos, const char *argName)
{
    char value[50];
    sprintf(value, "-%s", argName);

#ifdef DEBUG
    msgf("-%s ", argName);
#endif

    allocAndCopy(argv, pos, value);
}

void charArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    putArgName(argv, pos, argName);

    const char *_in = CHAR(asChar(in));
    allocAndCopy(argv, pos, _in);
#ifdef DEBUG
    msgf("%s\n", _in);
#endif
}

void logicalArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    if (asLogical(in))
    {
        putArgName(argv, pos, argName);
#ifdef DEBUG
        msgf("true\n");
#endif
    }
    else
    {
#ifdef DEBUG
        msgf("false\n");
#endif
    }
}

void integerArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    putArgName(argv, pos, argName);
    char value[50];
    sprintf(value, "%d", asInteger(in));

#ifdef DEBUG
    msgf("%d\n", asInteger(in));
#endif

    allocAndCopy(argv, pos, value);
}

void realArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    putArgName(argv, pos, argName);
    char value[50];
    sprintf(value, "%.14lf", asReal(in));

#ifdef DEBUG
    msgf("%.14lf\n", asReal(in));
#endif

    allocAndCopy(argv, pos, value);
}

void integerArrayArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    putArgName(argv, pos, argName);

    int size = LENGTH(in);

    for (int i = 0; i < size; i++)
    {
        char value[50];
        sprintf(value, "%d", INTEGER(in)[i]);
#ifdef DEBUG
        msgf("%d\n", INTEGER(in)[i]);
#endif
        allocAndCopy(argv, pos, value);
    }
}

void realArrayArg(char **argv, int *pos, const char *argName, SEXP in)
{
    if (isNull(in))
        return;

    putArgName(argv, pos, argName);

    int size = LENGTH(in);

    for (int i = 0; i < size; i++)
    {
        char value[50];
        sprintf(value, "%.14lf", REAL(in)[i]);

#ifdef DEBUG
        msgf("%.14lf\n", REAL(in)[i]);
#endif

        allocAndCopy(argv, pos, value);
    }
}