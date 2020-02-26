#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "gedisimulator/gediIO.h"
#include "tools/tools.h"
#include "argParse.h"

#define main addNoiseHDF
#define control addNoise_control
#define readCommands readCommands_addNoise
#define fprintf(stdout, ...) Rprintf(__VA_ARGS__)
    #include "gedisimulator/addNoiseHDF.c"
#undef fprintf
#undef readCommands
#undef control
#undef main


SEXP C_addNoiseHDF(
    // Input and output
    SEXP input,
    SEXP output,
    SEXP seed,
    SEXP linkNoise,
    SEXP dcBias,
    SEXP linkFsig,
    SEXP linkPsig,
    SEXP trueSig,
    SEXP bitRate)
{
    int argc = 1;
    char *argv[16];
    const char* algName = "addNoiseHDF";

    argv[0] = malloc((strlen(algName)+1) * sizeof(char));
    strcpy(argv[0], algName);

    PARSE_ARG(char, input);
    PARSE_ARG(char, output);
    PARSE_ARG(integer, seed);
    PARSE_ARG(realArray, linkNoise);
    PARSE_ARG(real, dcBias);
    PARSE_ARG(real, linkFsig);
    PARSE_ARG(real, linkPsig);
    PARSE_ARG(real, trueSig);
    PARSE_ARG(integer, bitRate);


#ifdef DEBUG   
    for (int i = 1; i < argc; i++) {
        Rprintf("%s ", argv[i]);
    }
#endif

    addNoiseHDF(argc, argv);

    for (int i = 0; i++; i < argc) {
        free(argv[i]);
    }

    
    return (ScalarInteger(0));
} 

#ifdef DEBUG
    #include "debug.c"

    int main() {
        initR();
        SEXP linkNoise = PROTECT(allocVector(REALSXP, 2));
        double *pLinkNoise = REAL(linkNoise);
        *(pLinkNoise++) = 3.0103;
        *(pLinkNoise++) = 0.95;
        pLinkNoise = NULL;

        C_addNoiseHDF(
            mkString("E:/Documentos/sample.h5"),
            mkString("E:/Documentos/sample_noised.h5"),
            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue);

            UNPROTECT(1);
    }
#endif
