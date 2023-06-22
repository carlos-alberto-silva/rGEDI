#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "gedisimulator/gediIO.h"
#include "tools/tools.h"
#include "argParse.h"

#define main collocateWaves_main
#define control collocateWaves_control
#define readCommands readCommands_collocateWaves
#include "gedisimulator/collocateWaves.c"
#undef readCommands
#undef control
#undef main

SEXP C_collocateWaves(
    // Input-output
    SEXP output,
    SEXP listAls,
    SEXP als,
    SEXP gedi,
    SEXP listGedi,
    SEXP readHDFgedi,
    SEXP lgw,
    SEXP readHDFlvis,
    SEXP bounds,
    SEXP leaveEmpty,
    SEXP lEPSG,
    SEXP aEPSG,
    SEXP beamList,
    SEXP skipBeams,
    SEXP readBeams,

    // Multi-mode options
    SEXP solveCofG,
    SEXP maxVshift,
    SEXP vStep,

    // Grid mode operation (multi-mode also applies)
    SEXP maxShift,
    SEXP step,

    // Optimiser mode operation (multi-mode also applies)
    SEXP simplex,
    SEXP anneal,
    SEXP fixFsig,
    SEXP geoError,
    SEXP quickGeo,
    SEXP optTol,
    SEXP maxIter,
    SEXP writeSimProg,
    SEXP writeWaves,
    SEXP nTriesAnneal,
    SEXP itersFixedT,
    SEXP kAnneal,
    SEXP tInitial,
    SEXP muAnneal,
    SEXP tMinAnneal,

    // Initial estimates. Will search around this point
    SEXP hOffset,
    SEXP offset,

    // Waveform characteristics
    SEXP fSigma,
    SEXP pSigma,
    SEXP readPulse,

    // Filters for input data
    SEXP minSense,
    SEXP maxZen,
    SEXP maxScanAng,
    SEXP minDense,
    SEXP decimate,
    SEXP noFilt,
    SEXP filtOut,
    SEXP smooth,
    SEXP checkCover,
    SEXP median,

    // Simulator settings. For simulator validation only
    SEXP noNorm,
    SEXP norm,
    SEXP allSimMeth,
    SEXP pulseBefore,

    // Octree to speed searching of ALS data. Not fully operational
    SEXP noOctree,
    SEXP octLevels,
    SEXP nOctPix)
{
    int argc = 1;
    char *argv[108];
    const char *algName = "collocateWaves";

    argv[0] = malloc((strlen(algName) + 1) * sizeof(char));
    strcpy(argv[0], algName);

    PARSE_ARG(char, output);
    PARSE_ARG(char, listAls);
    PARSE_ARG(char, als);
    PARSE_ARG(char, gedi);
    PARSE_ARG(char, listGedi);
    PARSE_ARG(logical, readHDFgedi);
    PARSE_ARG(logical, lgw);
    PARSE_ARG(logical, readHDFlvis);
    PARSE_ARG(realArray, bounds);
    PARSE_ARG(logical, leaveEmpty);
    PARSE_ARG(integer, lEPSG);
    PARSE_ARG(integer, aEPSG);
    PARSE_ARG(char, beamList);
    PARSE_ARG(char, skipBeams);
    PARSE_ARG(char, readBeams);

    PARSE_ARG(logical, solveCofG);
    PARSE_ARG(real, maxVshift);
    PARSE_ARG(real, vStep);

    PARSE_ARG(real, maxShift);
    PARSE_ARG(real, step);

    PARSE_ARG(logical, simplex);
    PARSE_ARG(logical, anneal);
    PARSE_ARG(logical, fixFsig);
    PARSE_ARG(realArray, geoError);
    PARSE_ARG(logical, quickGeo);
    PARSE_ARG(real, optTol);
    PARSE_ARG(integer, maxIter);
    PARSE_ARG(logical, writeSimProg);
    PARSE_ARG(char, writeWaves);
    PARSE_ARG(integer, nTriesAnneal);
    PARSE_ARG(integer, itersFixedT);
    PARSE_ARG(real, kAnneal);
    PARSE_ARG(real, tInitial);
    PARSE_ARG(real, muAnneal);
    PARSE_ARG(real, tMinAnneal);

    PARSE_ARG(realArray, hOffset);
    PARSE_ARG(real, offset);

    PARSE_ARG(real, fSigma);
    PARSE_ARG(real, pSigma);
    PARSE_ARG(char, readPulse);

    PARSE_ARG(real, minSense);
    PARSE_ARG(real, maxZen);
    PARSE_ARG(real, maxScanAng);
    PARSE_ARG(real, minDense);
    PARSE_ARG(real, decimate);
    PARSE_ARG(logical, noFilt);
    PARSE_ARG(real, filtOut);
    PARSE_ARG(real, smooth);
    PARSE_ARG(logical, checkCover);
    PARSE_ARG(logical, median);

    PARSE_ARG(logical, noNorm);
    PARSE_ARG(logical, norm);
    PARSE_ARG(logical, allSimMeth);
    PARSE_ARG(logical, pulseBefore);
    PARSE_ARG(logical, noOctree);
    PARSE_ARG(integer, octLevels);
    PARSE_ARG(integer, nOctPix);
    

#ifdef DEBUG
    Rprintf("./collocateWaves ");
    for (int i = 1; i < argc; i++)
    {
        Rprintf("%s ", argv[i]);
    }
    Rprintf("\n");
#endif

    return ScalarInteger(0);
}