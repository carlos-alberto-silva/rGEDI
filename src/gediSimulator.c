#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "gedisimulator/gediIO.h"
#include "tools/tools.h"
#include "argParse.h"

#define main gediRat
#define control metric_control
#define readCommands readCommands_rat
    #include "gedisimulator/gediRat.h"
    #include "gedisimulator/gediRat.c"
#undef readCommands
#undef control
#undef main



SEXP C_gediSimulator(
    // Input and output
    SEXP input,
    SEXP output,
    SEXP inList,
    SEXP ground,
    SEXP hdf,
    SEXP ascii,
    SEXP waveID,
    SEXP coord,
    SEXP listCoord,
    SEXP gridBound,
    SEXP gridStep,
    SEXP pSigma,
    SEXP pFWHM,
    SEXP readPulse,
    SEXP fSigma,
    SEXP wavefront,
    SEXP res,
    SEXP LVIS,
    SEXP topHat,
    SEXP sideLobe,
    SEXP lobeAng,
    SEXP checkCover,
    SEXP maxScanAng,
    SEXP decimate,
    SEXP pBuff,
    SEXP maxBins,
    SEXP countOnly,
    SEXP pulseAfter,
    SEXP pulseBefore,
    SEXP noNorm,
    SEXP noOctree,
    SEXP octLevels,
    SEXP nOctPix,
    SEXP decon,
    SEXP indDecon,
    SEXP readWave,
    SEXP listFiles,
    SEXP keepOld,
    SEXP useShadow,
    SEXP polyGround,
    SEXP nnGround,
    SEXP seed)
{
    int argc = 1;
    char *argv[108];
    const char* algName = "gediSimulator";

    argv[0] = malloc((strlen(algName)+1) * sizeof(char));
    strcpy(argv[0], algName);

    PARSE_ARG(char, input);
    PARSE_ARG(char, output);
    PARSE_ARG(char, inList);
    PARSE_ARG(logical, ground);
    PARSE_ARG(logical, hdf);
    PARSE_ARG(logical, ascii);
    PARSE_ARG(char, waveID);
    PARSE_ARG(realArray, coord);
    PARSE_ARG(char, listCoord);
    PARSE_ARG(realArray, gridBound);
    PARSE_ARG(real, gridStep);

    PARSE_ARG(real, pSigma);
    PARSE_ARG(real, pFWHM);
    PARSE_ARG(char, readPulse);
    PARSE_ARG(real, fSigma);
    PARSE_ARG(char, wavefront);
    PARSE_ARG(real, res);
    PARSE_ARG(logical, LVIS);
    PARSE_ARG(logical, topHat);
    PARSE_ARG(logical, sideLobe);
    PARSE_ARG(real, lobeAng);

    PARSE_ARG(logical, checkCover);
    PARSE_ARG(real, maxScanAng);
    PARSE_ARG(real, decimate);

    PARSE_ARG(real, pBuff);
    PARSE_ARG(integer, maxBins);
    PARSE_ARG(logical, countOnly);
    PARSE_ARG(logical, pulseAfter);
    PARSE_ARG(logical, pulseBefore);
    PARSE_ARG(logical, noNorm);

    PARSE_ARG(logical, noOctree);
    PARSE_ARG(integer, octLevels);
    PARSE_ARG(integer, nOctPix);

    PARSE_ARG(logical, decon);
    PARSE_ARG(logical, indDecon);
    PARSE_ARG(logical, readWave);

    PARSE_ARG(logical, listFiles);
    PARSE_ARG(logical, keepOld);
    PARSE_ARG(logical, useShadow);
    PARSE_ARG(logical, polyGround);
    PARSE_ARG(logical, nnGround);
    PARSE_ARG(integer, seed);


#ifdef DEBUG   
    for (int i = 1; i < argc; i++) {
        Rprintf("%s ", argv[i]);
    }
#endif

    gediRat(argc, argv);

    for (int i = 0; i++; i < argc) {
        free(argv[i]);
    }

    
    return (ScalarInteger(0));
} 

#ifdef DEBUG
    #include "debug.c"

    int main() {
        initR();
        SEXP gridBound = PROTECT(allocVector(REALSXP, 4));
        double *pGridBound = REAL(gridBound);
        *(pGridBound++) = 278204;
        *(pGridBound++) = 278296;
        *(pGridBound++) = 602203;
        *(pGridBound++) = 602296;
        pGridBound = NULL;

        C_gediSimulator(
            mkString("E:/Documentos/sample.las"),
            mkString("E:/Documentos/sample.h5"),
            R_NilValue,
            R_NilValue,
            ScalarLogical(1),
            R_NilValue, 
            R_NilValue,

            R_NilValue, 
            R_NilValue,
            gridBound,
            ScalarReal(15.0),

            ScalarReal(-1.0),
            ScalarReal(15.0),
            R_NilValue,
            ScalarReal(5.5),
            R_NilValue,
            ScalarReal(0.15),
            R_NilValue,
            R_NilValue,
            R_NilValue,
            ScalarReal(0.0),

            R_NilValue,
            ScalarReal(1000000.0),
            ScalarReal(1.0),

            ScalarInteger(20000000),
            ScalarInteger(1024),
            R_NilValue,
            R_NilValue,
            ScalarLogical(1),
            R_NilValue,

            R_NilValue,
            ScalarInteger(0),
            ScalarInteger(40),

            R_NilValue,
            R_NilValue,
            R_NilValue,

            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue,
            R_NilValue,
            ScalarInteger(200));

            UNPROTECT(1);
    }
#endif
