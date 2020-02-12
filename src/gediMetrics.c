#include <R.h>
#include <Rinternals.h>
#include "gedisimulator/gediIO.h"
#include "gedisimulator/gediMetric.h"
#include "tools/tools.h"

/*element reflectance*/
float rhoG;
float rhoC;

float *REALSXP2pFloat(SEXP);
char *INTSXP2pChar(SEXP);

char *INTSXP2pChar(SEXP x)
{
    if (isNull(x))
        return NULL;

    char *result;
    int size = LENGTH(x);
    result = (char *)malloc(sizeof(result) * size);
    for (int i = 0; i < size; i++)
    {
        result[i] = (char)INTEGER(x)[i];
    }
    return result;
};

metric_control *metric_makeControl(
    // Input and output
    const char *input,
    char writeFit,
    char writeGauss,
    char readBinLVIS,
    char readHDFlvis,
    char readHDFgedi,
    const char *level2,
    double *bounds,

    // Switches
    char ground,
    char useInt,
    char useFrac,
    float rhRes,
    float laiRes,
    float laiH,
    char noRHgauss,
    float gTol,
    float fhdHistRes,
    char forcePsigma,
    char bayesGround,
    char dontTrustGround,
    char noRoundCoord,

    // Adding noise
    float dcBias,
    float nSig,
    int *seed,
    float hNoise,
    char *linkNoise,
    float linkFsig,
    float linkPsig,
    float trueSig,
    int *bitRate,
    float maxDN,
    char renoise,
    float newPsig,
    float oldPsig,
    float addDrift,
    char missGround,
    float *minGap,

    // Photon counting
    char photonCount,
    float nPhotons,
    float photonWind,
    float noiseMult,

    // Denoising
    float meanN,
    float thresh,
    char varNoise,
    float varScale,
    float statsLen,
    char noiseTrack,
    float sWidth,
    float psWidth,
    float msWidth,
    char preMatchF,
    char postMatchF,
    const char *pFile,
    float gWidth,
    float minGsig,
    float minWidth,
    char medNoise,
    char varDrift,
    float *driftFac,
    float varRhoG,
    float varRhoC,
    float pSigma,
    char gold,
    float deconTol)
{
    metric_control *dimage = NULL;

    /*allocate structures*/
    if (!(dimage = (metric_control *)calloc(1, sizeof(metric_control))))
    {
        fprintf(stderr, "error metric_control allocation.\n");
        exit(1);
    }
    if (!(dimage->gediIO.den = (denPar *)calloc(1, sizeof(denPar))))
    {
        fprintf(stderr, "error metric_control allocation.\n");
        exit(1);
    }
    if (!(dimage->gediIO.gFit = (denPar *)calloc(1, sizeof(denPar))))
    {
        fprintf(stderr, "error metric_control allocation.\n");
        exit(1);
    }

    /*defaults*/
    /*input/output*/
    dimage->gediIO.nFiles = 1;
    dimage->gediIO.inList = chChalloc(dimage->gediIO.nFiles, "inList", 0);
    dimage->gediIO.inList[0] = challoc(200, "inList", 0);
    strcpy(&(dimage->gediIO.inList[0][0]), "/Users/stevenhancock/data/gedi/analysis/side_lobe/laselva/contLaselva.0.12.wave");
    strcpy(dimage->outRoot, "teastMetric");
    dimage->maxGauss = 20;
    dimage->opooMet = NULL;
    dimage->opooGauss = NULL;
    dimage->hdfGedi = NULL;
    /*scan settings*/
    dimage->gediIO.pSigma = 0.764331; /*pulse length*/
    dimage->gediIO.fSigma = 5.5;      /*footprint width*/
    dimage->gediIO.res = 0.15;

    /*switches*/
    dimage->writeFit = 0;
    dimage->gediIO.ground = 0;
    dimage->gediIO.useInt = 0;
    dimage->gediIO.useCount = 1;
    dimage->gediIO.useFrac = 0;
    dimage->rhRes = 5.0;
    dimage->laiRes = 10.0;
    dimage->maxLAIh = 30.0;
    dimage->bayesGround = 0;
    dimage->noise.missGround = 0;
    dimage->noise.linkNoise = 0;
    dimage->noise.driftFact = 0.0;
    dimage->gediIO.linkPsig = 0.764331; /*pulse length*/
    dimage->gediIO.linkFsig = 5.5;      /*footprint width*/
    dimage->noise.trueSig = 5.0;
    dimage->noise.deSig = 0.0; //0.1; //4.0*0.15/2.355;
    dimage->noise.bitRate = 12;
    dimage->noise.maxDN = 4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
    dimage->noise.minGap = 0.0;
    dimage->noRHgauss = 0;              /*do find RH metrics by Gaussian fitting*/
    dimage->renoiseWave = 0;            /*do not denoise "truth"*/
    dimage->noise.newPsig = -1.0;       /*leave blank*/
    dimage->gediIO.dontTrustGround = 0; /*do trust ground in waveforms, if there*/
    dimage->readBinLVIS = 0;            /*read ASCII rather than binary LVIS*/
    dimage->readHDFlvis = 0;            /*read ASCII rather than HDF5 LVIS*/
    dimage->readHDFgedi = 0;            /*read ASCII rather than HDF5 GEDI*/
    dimage->gediIO.readPsigma = 1;      /*read pSigma from file*/
    dimage->coord2dp = 1;               /*round up coords in output*/
    dimage->useBounds = 0;              /*process all data provided*/
    dimage->writeGauss = 0;             /*do not write Gaussian parameters*/

    /*set default denoising parameters*/
    setDenoiseDefault(dimage->gediIO.den);
    dimage->gediIO.den->meanN = 0.0;         /*we haven't added noise yet*/
    dimage->gediIO.den->thresh = 0.00000001; /*tiny number as no noise yet*/
    dimage->gediIO.den->noiseTrack = 0;
    dimage->gediIO.den->minWidth = 0;
    dimage->gediIO.den->varNoise = 0;
    dimage->gediIO.den->threshScale = 1.5;
    dimage->gediIO.den->fitGauss = 0;
    dimage->gediIO.den->psWidth = 0.0;

    /*set default Gaussian fitting  parameters*/
    setDenoiseDefault(dimage->gediIO.gFit);
    dimage->gediIO.gFit->meanN = 0.0;          /*no denoising here*/
    dimage->gediIO.gFit->thresh = 0.000000005; /*no denoising here*/
    dimage->gediIO.gFit->noiseTrack = 0;       /*no denoising here*/
    dimage->gediIO.gFit->minWidth = 0;         /*no denoising here*/
    dimage->gediIO.gFit->varNoise = 0;         /*no denoising here*/
    dimage->gediIO.gFit->gWidth = 1.2;
    dimage->gediIO.gFit->sWidth = 0.0;
    dimage->gediIO.gFit->fitGauss = 1;
    dimage->gediIO.gFit->minGsig = 0.764331;
    /*noise parameters*/
    dimage->noise.meanN = 0.0;
    dimage->noise.nSig = nSig;
    dimage->bThresh = 0.001;
    dimage->noise.hNoise = hNoise;
    dimage->noise.offset = 94.0;
    /*projection, not yet used*/
    dimage->gediIO.wEPSG = 4326; /*waveforms*/
    dimage->gediIO.bEPSG = 4326; /*bounds*/
    /*LVIS data*/
    dimage->lvis.data = NULL;
    dimage->hdfLvis = NULL;
    /*LVIS level2 data*/
    dimage->readL2 = 0; /*do not read L2*/
    /*photon counting*/
    dimage->ice2 = 0; /*GEDI mode, rather than ICESat-2*/

    // Change parameters based on user input

    // TODO if readHDFlvis
    dimage->readHDFlvis = 1;
    dimage->gediIO.readPsigma = 0;

    // TODO if varNoise
    dimage->gediIO.den->varNoise = 1;

    // TODO if varScale
    dimage->gediIO.den->varNoise = 1;
    dimage->gediIO.den->threshScale = 3;

    // TODO minWidth
    dimage->gediIO.den->minWidth = 2;

    // TODO sWidth
    dimage->gediIO.den->sWidth = 8;

    // TODO statsLen
    dimage->gediIO.den->statsLen = 10;

#ifdef USEPHOTON
    dimage->photonCount.designval = 2.1;
    dimage->photonCount.prob = NULL;
    dimage->photonCount.pBins = 0;
    dimage->photonCount.H = 200.0;
    dimage->photonCount.noise_mult = 0.1;
    dimage->photonCount.rhoVrhoG = 1.0;
    dimage->photonCount.writeHDF = 0; /*write ASCII by default*/
    dimage->photonCount.hdf = NULL;
#endif
    /*others*/
    rhoG = 0.4;
    rhoC = 0.57;
    dimage->rhoRatio = rhoC / rhoG;
    dimage->gTol = 0.0;
    dimage->gediIO.nMessages = 200;
    dimage->fhdHistRes = 0.001;

    /*read the command line*/
    TTIDY((void **)dimage->gediIO.inList, dimage->gediIO.nFiles);
    dimage->gediIO.inList = NULL;
    dimage->gediIO.nFiles = 1;
    dimage->gediIO.inList = chChalloc(dimage->gediIO.nFiles, "input name list", 0);
    dimage->gediIO.inList[0] = challoc((uint64_t)strlen(input) + 1, "input name list", 0);
    strcpy(dimage->gediIO.inList[0], input);
    if (level2)
    {
        dimage->readL2 = 1;
        strcpy(dimage->l2namen, level2);
    }
    if (writeFit)
        dimage->writeFit = 1;
    if (writeGauss)
        dimage->writeGauss = 1;
    if (dcBias)
        dimage->noise.meanN = dimage->noise.offset = dcBias;
    if (ground)
        dimage->gediIO.ground = 1;
    if (varNoise)
        dimage->gediIO.den->varNoise = 1;
    if (medNoise)
        dimage->gediIO.den->medStats = 1;
    if (noiseTrack)
        dimage->gediIO.den->noiseTrack = 1;
    if (pFile)
    {
        strcpy(dimage->gediIO.den->pNamen, pFile);
        dimage->gediIO.den->deconGauss = 0;
    }
    if (gold)
        dimage->gediIO.den->deconMeth = 0;
    if (preMatchF)
        dimage->gediIO.den->preMatchF = 1;
    if (postMatchF)
        dimage->gediIO.den->posMatchF = 1;
    if (useInt)
    {
        dimage->gediIO.useInt = 1;
        dimage->gediIO.useCount = 0;
        dimage->gediIO.useFrac = 0;
    }
    if (useFrac)
    {
        dimage->gediIO.useInt = 0;
        dimage->gediIO.useCount = 0;
        dimage->gediIO.useFrac = 1;
    }
    if (linkNoise)
    {
        dimage->noise.linkNoise = 1;
        dimage->noise.linkM = linkNoise[0];
        dimage->noise.linkCov = linkNoise[1];
    }
    if (missGround)
        dimage->noise.missGround = 1;
    if (minGap)
    {
        dimage->noise.missGround = 1;
        dimage->noise.minGap = *minGap;
    }
    if (bayesGround)
        dimage->bayesGround = 1;
    if (bitRate)
    {
        dimage->noise.bitRate = *bitRate;
        dimage->noise.maxDN = pow(2.0, (float)dimage->noise.bitRate);
    }
    if (noRHgauss)
        dimage->noRHgauss = 1;
    if (renoise)
        dimage->renoiseWave = 1;
    if (dontTrustGround)
        dimage->gediIO.dontTrustGround = 1;
    if (readBinLVIS)
        dimage->readBinLVIS = 1;
    if (readHDFlvis)
    {
        dimage->readHDFlvis = 1;
        dimage->gediIO.readPsigma = 0;
    }
    if (readHDFgedi)
        dimage->readHDFgedi = 1;
    if (forcePsigma)
        dimage->gediIO.readPsigma = 0;
    if (noRoundCoord)
        dimage->coord2dp = 0;
    if (bounds)
    {
        dimage->useBounds = 1;
        dimage->minX = dimage->gediIO.bounds[0] = bounds[0];
        dimage->minY = dimage->gediIO.bounds[1] = bounds[1];
        dimage->maxX = dimage->gediIO.bounds[2] = bounds[2];
        dimage->maxY = dimage->gediIO.bounds[3] = bounds[3];
    }
    if (varDrift)
    {
        dimage->gediIO.den->corrDrift = 1;
        dimage->gediIO.den->varDrift = 1;
    }
    if (driftFac)
    {
        dimage->gediIO.den->corrDrift = 1;
        dimage->gediIO.den->varDrift = 0;
        dimage->gediIO.den->fixedDrift = *driftFac;
    }
    if (seed)
        srand(*seed);

#ifdef USEPHOTON
    if (photonCount)
    {
        dimage->ice2 = 1;
    }
    if (nPhotons)
    {
        dimage->photonCount.designval = nPhotons;
    }
    if (photonWind)
    {
        dimage->photonCount.H = photonWind;
    }
    if (noiseMult)
    {
        dimage->photonCount.noise_mult = noiseMult;
    }
#endif

    /*read deconvolution pulse if needed*/
    if (dimage->gediIO.den->preMatchF || dimage->gediIO.den->preMatchF || dimage->gediIO.den->deconMeth >= 0)
        readPulse(dimage->gediIO.den);
    if ((!dimage->gediIO.ground) && (dimage->noise.missGround))
    {
        fprintf(stderr, "Noise option conflict. Cannot use missGround without ground\n");
        exit(1);
    }

    return dimage;
}

float *REALSXP2pFloat(SEXP x)
{
    if (isNull(x))
        return NULL;

    float *result;
    int size = LENGTH(x);
    result = (float *)malloc(sizeof(result) * size);
    for (int i = 0; i < size; i++)
    {
        result[i] = (float)REAL(x)[i];
    }
    return result;
}

#define exit(n) return (ScalarInteger(n))
SEXP C_gediMetrics(
    // Input and output
    SEXP input,
    SEXP writeFit,
    SEXP writeGauss,
    SEXP readBinLVIS,
    SEXP readHDFlvis,
    SEXP readHDFgedi,
    SEXP level2,
    SEXP bounds,

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

    // Photon counting
    SEXP photonCount,
    SEXP nPhotons,
    SEXP photonWind,
    SEXP noiseMult,

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
    SEXP gWidth,
    SEXP minGsig,
    SEXP minWidth,
    SEXP medNoise,
    SEXP varDrift,
    SEXP driftFac,
    SEXP varRhoG,
    SEXP varRhoC,
    SEXP pSigma,
    SEXP gold,
    SEXP deconTol)
{
    int i = 0;
    metric_control *dimage = NULL;
    dataStruct *data = NULL;
    metStruct *metric = NULL;
    void setL2ground(dataStruct *, int, metric_control *);
    void findMetrics(metStruct *, float *, int, float *, float *, int, double *, metric_control *, dataStruct *);
    void tidySMoothPulse();
    void alignElevation(double, double, float *, int);
    void writeResults(dataStruct *, metric_control *, metStruct *, int, float *, float *, char *);
    void determineTruth(dataStruct *, metric_control *);
    void modifyTruth(dataStruct *, noisePar *);
    void checkWaveformBounds(dataStruct *, metric_control *);
    void photonCountCloud(float *, dataStruct *, photonStruct *, char *, int, denPar *, noisePar *);
    float *processed = NULL, *denoised = NULL;
    ;
    int j = 0;

    /*read command Line*/
    dimage = metric_makeControl(
        // Input and output
        CHAR(asChar(input)),
        (char)asLogical(writeFit),
        (char)asLogical(writeGauss),
        (char)asLogical(readBinLVIS),
        (char)asLogical(readHDFlvis),
        (char)asLogical(readHDFgedi),
        isNull(level2) ? NULL : CHAR(asChar(level2)),
        isNull(bounds) ? NULL : REAL(bounds),

        // Switches
        (char)asLogical(ground),
        (char)asLogical(useInt),
        (char)asLogical(useFrac),
        (float)asReal(rhRes),
        (float)asReal(laiRes),
        (float)asReal(laiH),
        (char)asLogical(noRHgauss),
        (float)asReal(gTol),
        (float)asReal(fhdHistRes),
        (char)asLogical(forcePsigma),
        (char)asLogical(bayesGround),
        (char)asLogical(dontTrustGround),
        (char)asLogical(noRoundCoord),

        // Adding noise
        (float)asReal(dcBias),
        (float)asReal(nSig),
        isNull(seed) ? NULL : INTEGER(seed),
        (float)asReal(hNoise),
        INTSXP2pChar(linkNoise),
        (float)asReal(linkFsig),
        (float)asReal(linkPsig),
        (float)asReal(trueSig),
        isNull(bitRate) ? NULL : INTEGER(bitRate),
        (float)asReal(maxDN),
        (char)asLogical(renoise),
        (float)asReal(newPsig),
        (float)asReal(oldPsig),
        (float)asReal(addDrift),
        (char)asLogical(missGround),
        REALSXP2pFloat(minGap),

        // Photon counting
        (char)asLogical(photonCount),
        (float)asReal(nPhotons),
        (float)asReal(photonWind),
        (float)asReal(noiseMult),

        // Denoising
        (float)asReal(meanN),
        (float)asReal(thresh),
        (char)asLogical(varNoise),
        (float)asReal(varScale),
        (float)asReal(statsLen),
        (char)asLogical(noiseTrack),
        (float)asReal(sWidth),
        (float)asReal(psWidth),
        (float)asReal(msWidth),
        (char)asLogical(preMatchF),
        (char)asLogical(postMatchF),
        isNull(pFile) ? NULL : CHAR(asChar(pFile)),
        (float)asReal(gWidth),
        (float)asReal(minGsig),
        (float)asReal(minWidth),
        (char)asLogical(medNoise),
        (char)asLogical(varDrift),
        REALSXP2pFloat(driftFac),
        (float)asReal(varRhoG),
        (float)asReal(varRhoC),
        (float)asReal(pSigma),
        (char)asLogical(gold),
        (float)asReal(deconTol));

    /*set link noise if needed*/
    dimage->noise.linkSig = setNoiseSigma(
        dimage->noise.linkM,
        dimage->noise.linkCov,
        dimage->gediIO.linkPsig,
        dimage->gediIO.linkFsig,
        rhoC,
        rhoG);

/*set photon rates if needed*/
#ifdef USEPHOTON
    if (dimage->ice2 || dimage->pclPhoton)
        setPhotonRates(&dimage->photonCount);
#endif

    /*allocate metric array*/
    if (!(metric = (metStruct *)calloc(1, sizeof(metStruct))))
    {
        fprintf(stderr, "error metric structure allocation.\n");
        exit(1);
    }

    /*loop over files*/
    for (i = 0; i < dimage->gediIO.nFiles; i++)
    {
        if ((i % dimage->gediIO.nMessages) == 0)
            fprintf(
                stdout,
                "Wave %d of %d\n",
                i + 1,
                dimage->gediIO.nFiles);

        /*read waveform*/
        if (dimage->readBinLVIS)
            data = readBinaryLVIS(dimage->gediIO.inList[0],
                                  &dimage->lvis,
                                  i,
                                  &dimage->gediIO);
        else if (dimage->readHDFlvis)
            data = unpackHDFlvis(dimage->gediIO.inList[0],
                                 &dimage->hdfLvis,
                                 &dimage->gediIO,
                                 i);
        else if (dimage->readHDFgedi)
            data = unpackHDFgedi(dimage->gediIO.inList[0],
                                 &dimage->gediIO,
                                 &dimage->hdfGedi,
                                 i);
        else
            data = readASCIIdata(dimage->gediIO.inList[i],
                                 &(dimage->gediIO));
        if (dimage->readL2)
            setL2ground(data, i, dimage);

        /*check bounds if needed*/
        if (dimage->useBounds)
            checkWaveformBounds(data, dimage);

        /*is the data usable*/
        if (data->usable)
        {
            /*denoise and change pulse if needed*/
            if (dimage->renoiseWave)
                modifyTruth(data, &dimage->noise);

            /*determine truths before noising*/
            determineTruth(data, dimage);

            /*add noise if needed*/
            if (!dimage->pclPhoton)
                addNoise(data, &dimage->noise,
                         dimage->gediIO.fSigma,
                         dimage->gediIO.pSigma,
                         dimage->gediIO.res,
                         rhoC,
                         rhoG);
            else
                data->noised = uncompressPhotons(
                    data->wave[data->useType],
                    data,
                    &dimage->photonCount,
                    &dimage->noise,
                    &dimage->gediIO);

            /*process waveform*/
            /*denoise, or*if we are doing PCL on photon counting, convert to photon count*/
            denoised = processFloWave(data->noised,
                                      data->nBins,
                                      dimage->gediIO.den,
                                      1.0);
            if (dimage->pclPhoton)
                for (j = 0; j < data->nBins; j++)
                    fprintf(stdout,
                            "wave %d %d %f %f %f\n",
                            i,
                            j,
                            denoised[j],
                            data->wave[0][j],
                            data->noised[j]);

            /*check that the wave is still usable*/
            if (checkUsable(denoised, data->nBins))
            {
                /*are we in GEDI mode?*/
                if (!dimage->ice2)
                {

                    /*Gaussian fit*/
                    if (dimage->noRHgauss == 0)
                        processed = processFloWave(
                            denoised,
                            data->nBins,
                            dimage->gediIO.gFit,
                            1.0);

                    /*shift Gaussian centres to align to absolute elevation*/
                    alignElevation(data->z[0],
                                   data->z[data->nBins - 1],
                                   dimage->gediIO.gFit->gPar,
                                   dimage->gediIO.gFit->nGauss);

                    /*determine metrics*/
                    findMetrics(metric,
                                dimage->gediIO.gFit->gPar,
                                dimage->gediIO.gFit->nGauss,
                                denoised,
                                data->noised,
                                data->nBins,
                                data->z,
                                dimage,
                                data);

                    /*write results*/
                    if (dimage->readBinLVIS || dimage->readHDFlvis || dimage->readHDFgedi)
                        writeResults(data,
                                     dimage,
                                     metric,
                                     i,
                                     denoised,
                                     processed,
                                     dimage->gediIO.inList[0]);
                    else
                        writeResults(data,
                                     dimage,
                                     metric,
                                     i,
                                     denoised,
                                     processed,
                                     dimage->gediIO.inList[i]);
                }
                else
                { /*ICESat-2 mode*/
                    photonCountCloud(denoised,
                                     data,
                                     &dimage->photonCount,
                                     dimage->outRoot,
                                     i,
                                     dimage->gediIO.den,
                                     &dimage->noise);
                } /*operation mode switch*/
            }
            else
            { /*still usable after denoising?*/
                fprintf(stderr, "No longer usable\n");
            }
        } /*is the data usable*/

        /*tidy as we go along*/
        TIDY(processed);
        TIDY(denoised);
        if (data)
        {
            TIDY(data->noised);
            if (dimage->readHDFgedi)
            { /*pointer to array. do not free*/
                data->wave[0] = NULL;
                if (data->ground)
                    data->ground[0] = NULL;
            }
            TTIDY((void **)data->ground, data->nWaveTypes);
            TTIDY((void **)data->wave, data->nWaveTypes);
            TIDY(data->totE);
            TIDY(data->z);
            TIDY(data);
        }
        TIDY(dimage->gediIO.gFit->gPar);
        TIDY(dimage->gediIO.den->gPar);
        dimage->gediIO.den->nGauss = 0;
        dimage->gediIO.gFit->nGauss = 0;
        TIDY(metric->rhMax);
        TIDY(metric->rhInfl);
        TIDY(metric->rhReal);
        TIDY(metric->rh);
        TIDY(metric->bGr);
        TIDY(metric->tLAI);
        TIDY(metric->gLAI);
        TIDY(metric->hgLAI);
        TIDY(metric->hiLAI);
        TIDY(metric->hmLAI);
        //TIDY(metric->LmomGau);
        //TIDY(metric->LmomRea);
        //TIDY(metric->LmomInf);
        //TIDY(metric->LmomMax);
    } /*file loop*/

    /*TIDY LVIS data if it was read*/
    if (dimage->readBinLVIS)
        TIDY(dimage->lvis.data);
    if (dimage->readHDFgedi)
        dimage->hdfGedi = tidyGediHDF(dimage->hdfGedi);

    if (dimage->writeGauss)
        fprintf(stdout,
                "Written to %s.gauss.txt\n",
                dimage->outRoot);
    if (!dimage->ice2)
        fprintf(stdout,
                "Written to %s.metric.txt\n",
                dimage->outRoot);
#ifdef USEPHOTON
    else
        fprintf(stdout,
                "Written to %s\n",
                dimage->photonCount.outNamen);
#endif

    /*tidy up arrays*/
    tidySMoothPulse();
    TIDY(metric);
    if (dimage)
    {
        if (dimage->lvisL2)
        {
            TIDY(dimage->lvisL2->lfid);
            TIDY(dimage->lvisL2->shotN);
            TIDY(dimage->lvisL2->zG);
            TIDY(dimage->lvisL2);
        }
        if (dimage->readBinLVIS || dimage->readHDFlvis || dimage->readHDFgedi)
            TTIDY((void **)dimage->gediIO.inList, 1);
        else
            TTIDY((void **)dimage->gediIO.inList,
                  dimage->gediIO.nFiles);
        dimage->gediIO.inList = NULL;
        TIDY(dimage->gediIO.noiseSigs.threshN);
        TIDY(dimage->gediIO.noiseSigs.threshS);
        TIDY(dimage->gediIO.noiseSigs.probNoise);
        TIDY(dimage->gediIO.noiseSigs.probMiss);
        if (dimage->opooMet)
        {
            fclose(dimage->opooMet);
            dimage->opooMet = NULL;
        }
        if (dimage->opooGauss)
        {
            fclose(dimage->opooGauss);
            dimage->opooGauss = NULL;
        }
#ifdef USEPHOTON
        if (dimage->photonCount.opoo)
        {
            fclose(dimage->photonCount.opoo);
            dimage->photonCount.opoo = NULL;
        }
        TIDY(dimage->photonCount.prob);
#endif
        if (dimage->gediIO.den)
        {
            TTIDY((void **)dimage->gediIO.den->pulse, 2);
            TIDY(dimage->gediIO.den->matchPulse);
            TIDY(dimage->gediIO.den->hardPulse);
            TIDY(dimage->gediIO.den);
        }
        if (dimage->gediIO.gFit)
        {
            TTIDY((void **)dimage->gediIO.gFit->pulse, 2);
            TIDY(dimage->gediIO.gFit);
        }
        dimage->hdfLvis = tidyLVISstruct(dimage->hdfLvis);
        TIDY(dimage);
    }
    return (ScalarInteger(0));
} /*main*/

#undef exit