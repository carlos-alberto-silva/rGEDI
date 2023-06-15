#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "gedisimulator/gediIO.h"
#include "tools/tools.h"
#include "argParse.h"

#define main gediMetric
#define control metric_control
#define readCommands readCommands_metric
    #include "gedisimulator/gediMetric.h"
    #include "gedisimulator/gediMetric.c"
#undef readCommands
#undef control
#undef main


SEXP C_gediMetrics(
    // Input output
    SEXP input,
    SEXP outRoot,
    SEXP inList,
    SEXP writeFit,
    SEXP writeGauss,
    SEXP readBinLVIS,
    SEXP readHDFlvis,
    SEXP readHDFgedi,
    SEXP level2,
    SEXP bounds,
    SEXP beamList,
    SEXP skipBeams,
    SEXP readBeams,

    // Switches
    SEXP ground,
    SEXP useInt,
    SEXP useFrac,
    SEXP rhRes,
    SEXP laiRes,
    SEXP laiH,
    SEXP noRHgauss,
    SEXP gTol,
    SEXP fhdHistRes,
    SEXP forcePsigma,
    SEXP bayesGround,
    SEXP dontTrustGround,
    SEXP noRoundCoord,
    SEXP noCanopy,

    // Adding noise
    SEXP dcBias,
    SEXP nSig,
    SEXP seed,
    SEXP hNoise,
    SEXP linkNoise,
    SEXP linkFsig,
    SEXP linkPsig,
    SEXP trueSig,
    SEXP bitRate,
    SEXP maxDN,
    SEXP renoise,
    SEXP newPsig,
    SEXP oldPsig,
    SEXP addDrift,
    SEXP missGround,
    SEXP minGap,

    // Photon Counting
    SEXP photonCount,
    SEXP pcl,
    SEXP nPhotons,
    SEXP photonWind,
    SEXP noiseMult,
    SEXP rhoVrhoG,
    SEXP nPhotC,
    SEXP nPhotG,
    SEXP photHDF,

    // Denoising
    SEXP meanN,
    SEXP thresh,
    SEXP varNoise,
    SEXP varScale,
    SEXP statsLen,
    SEXP noiseTrack,
    SEXP sWidth,
    SEXP psWidth,
    SEXP msWidth,
    SEXP preMatchF,
    SEXP postMatchF,
    SEXP pFile,
    SEXP aggParams)
{
    int argc = 1;
    char *argv[108];
    const char* algName = "gediMetric";

    argv[0] = malloc((strlen(algName)+1) * sizeof(char));
    strcpy(argv[0], algName);

    int i = 0;
    #define paramsFromList VECTOR_ELT(aggParams, i++)
    SEXP gWidth = paramsFromList;
    SEXP minGsig = paramsFromList;
    SEXP minWidth = paramsFromList;
    SEXP medNoise = paramsFromList;
    SEXP varDrift = paramsFromList;
    SEXP driftFac = paramsFromList;
    SEXP rhoG = paramsFromList;
    SEXP rhoC = paramsFromList;
    SEXP pSigma = paramsFromList;
    SEXP gold = paramsFromList;
    SEXP deconTol = paramsFromList;
    #undef paramsFromList

    PARSE_ARG(char, input);
    PARSE_ARG(char, outRoot);
    PARSE_ARG(char, inList);
    PARSE_ARG(logical, writeFit);
    PARSE_ARG(logical, writeGauss);
    PARSE_ARG(logical, readBinLVIS);
    PARSE_ARG(logical, readHDFlvis);
    PARSE_ARG(logical, readHDFgedi);
    PARSE_ARG(char, level2);
    PARSE_ARG(realArray, bounds);
    PARSE_ARG(char, beamList);
    PARSE_ARG(char, skipBeams);
    PARSE_ARG(char, readBeams);
    PARSE_ARG(logical, ground);
    PARSE_ARG(logical, useInt);
    PARSE_ARG(logical, useFrac);
    PARSE_ARG(real, rhRes);
    PARSE_ARG(real, laiRes);
    PARSE_ARG(real, laiH);
    PARSE_ARG(logical, noRHgauss);
    PARSE_ARG(real, gTol);
    PARSE_ARG(real, fhdHistRes);
    PARSE_ARG(logical, forcePsigma);
    PARSE_ARG(logical, bayesGround);
    PARSE_ARG(logical, dontTrustGround);
    PARSE_ARG(logical, noRoundCoord);
    PARSE_ARG(logical, noCanopy);
    PARSE_ARG(real, dcBias);
    PARSE_ARG(real, nSig);
    PARSE_ARG(integer, seed);
    PARSE_ARG(real, hNoise);
    PARSE_ARG(realArray, linkNoise);
    PARSE_ARG(real, linkFsig);
    PARSE_ARG(real, linkPsig);
    PARSE_ARG(real, trueSig);
    PARSE_ARG(integer, bitRate);
    PARSE_ARG(real, maxDN);
    PARSE_ARG(logical, renoise);
    PARSE_ARG(real, newPsig);
    PARSE_ARG(real, oldPsig);
    PARSE_ARG(real, addDrift);
    PARSE_ARG(logical, missGround);
    PARSE_ARG(real, minGap);
    PARSE_ARG(logical, photonCount);
    PARSE_ARG(logical, pcl);
    PARSE_ARG(real, nPhotons);
    PARSE_ARG(real, photonWind);
    PARSE_ARG(real, noiseMult);
    PARSE_ARG(real, rhoVrhoG);
    PARSE_ARG(real, nPhotC);
    PARSE_ARG(real, nPhotG);
    PARSE_ARG(logical, photHDF);
    PARSE_ARG(real, meanN);
    PARSE_ARG(real, thresh);
    PARSE_ARG(logical, varNoise);
    PARSE_ARG(real, varScale);
    PARSE_ARG(real, statsLen);
    PARSE_ARG(logical, noiseTrack);
    PARSE_ARG(real, sWidth);
    PARSE_ARG(real, psWidth);
    PARSE_ARG(real, msWidth);
    PARSE_ARG(logical, preMatchF);
    PARSE_ARG(logical, postMatchF);
    PARSE_ARG(char, pFile);
    PARSE_ARG(real, gWidth);
    PARSE_ARG(real, minGsig);
    PARSE_ARG(real, minWidth);
    PARSE_ARG(logical, medNoise);
    PARSE_ARG(logical, varDrift);
    PARSE_ARG(real, driftFac);
    PARSE_ARG(real, rhoG);
    PARSE_ARG(real, rhoC);
    PARSE_ARG(real, pSigma);
    PARSE_ARG(logical, gold);
    PARSE_ARG(real, deconTol);


#ifdef DEBUG   
    Rprintf("./gediMetric ");
    for (int i = 1; i < argc; i++) {
        Rprintf("%s ", argv[i]);
    }
    Rprintf("\n");
#endif

    GetRNGstate();
    int status = gediMetric(argc, argv);
    if (status != 0) {
        REprintf("Something went wrong!");
    }
    PutRNGstate();

    for (int i = 0; i < argc; i++) {
        free(argv[i]);
    }

    
    return (ScalarInteger(0));
} 