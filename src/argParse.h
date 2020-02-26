#ifndef RGEDI_ARGPARSE_H
#define RGEDI_ARGPARSE_H
#include <Rinternals.h>

#define PARSE_ARG(typeName, n) typeName##Arg(argv, &argc, #n, n)

void charArg(char **argv, int *pos, const char *argName, SEXP in);
void integerArg(char **argv, int *pos, const char *argName, SEXP in);
void integerArrayArg(char **argv, int *pos, const char *argName, SEXP in);
void logicalArg(char **argv, int *pos, const char *argName, SEXP in);
void putArgName(char **argv, int *pos, const char *argName);
void realArg(char **argv, int *pos, const char *argName, SEXP in);
void realArrayArg(char **argv, int *pos, const char *argName, SEXP in);
#endif /* RGEDI_ARGPARSE_H */
