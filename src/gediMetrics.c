#include <R.h>
#include <Rinternals.h>
#include "gedisimulator/gediMetric.h"

/*element reflectance*/
float rhoG;
float rhoC;

metric_control *metric_makeControl(
    const char* input,
    char writeFit,
    char writeGauss,
    char readBinLVIS,
    char readHDFlvis,
    char readHDFgedi,
    const char* level2,
    double* bounds,
    char ground,
    char useInt,
    char useFrac,
    float rhRes = 5,
    float laiRes = 10.0,
    float laiH = 30.0,
    char noRHgauss,
    float gTol = 0.0,
    float fhdHistRes = 0.001,
    char forcePsigma,
    char bayesGround,
    char dontTrustGround,
    char noRoundCoord,
    NumericVector dcBias = NumericVector(0),
    float nSig = 0.0,
    int seed = 1,
    float hNoise = 0.0,
    NumericVector linkNoise = NumericVector(0),
    float linkFsig = 5.5,
    float linkPsig = 0.764331,
    float trueSig = 5.0,
    int bitRate = 12,
    float maxDN = 4096.0,
    char renoise = false,
    float newPsig = 1.0,
    float oldPsig = 0.764331,
    float addDrift = 0.0,
    char missGround = 0,
    NumericVector minGap = NumericVector(0),
    char photonCount = false,
    float nPhotons = 2.1,
    float photonWind = 200,
    float noiseMult = 0.1,
    float meanN = 0.0,
    float thresh = 0.00000001,
    char varNoise = false,
    float varScale = 1.5,
    float statsLen = 30.0,
    char noiseTrack = false,
    float sWidth = 0.0,
    float psWidth = 0.0,
    float msWidth = 0.0,
    char preMatchF = false,
    char postMatchF = false,
    CharacterVector pFile = CharacterVector(0),
    float gWidth = 1.2,
    float minGsig = 0.764331,
    float minWidth = 0,
    char medNoise = false,
    char varDrift = false,
    NumericVector driftFac = NumericVector(0),
    float varRhoG = 0.4,
    float varRhoC = 0.57,
    float pSigma = 1.0,
    char gold = false,
    float deconTol = 0.0000001)
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
        dimage->noise.meanN = dimage->noise.offset = dcBias[0];
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
        strcpy(dimage->gediIO.den->pNamen, pFile[0]);
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
        dimage->noise.minGap = minGap[0];
    }
    if (bayesGround)
        dimage->bayesGround = 1;
    if (bitRate)
    {
        dimage->noise.bitRate = bitRate;
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
        dimage->gediIO.den->fixedDrift = driftFac[0];
#ifdef USEPHOTON
    }
    if (!strncasecmp(argv[i], "-photonCount", 12))
    {
        dimage->ice2 = 1;
    }
    if (!strncasecmp(argv[i], "-nPhotons", 9))
    {
        checkArguments(1, i, argc, "-nPhotons");
        dimage->photonCount.designval = atof(argv[++i]);
    }
    if (!strncasecmp(argv[i], "-photonWind", 11))
    {
        checkArguments(1, i, argc, "-photonWind");
        dimage->photonCount.H = atof(argv[++i]);
    }
    if (!strncasecmp(argv[i], "-noiseMult", 10))
    {
        checkArguments(1, i, argc, "-noiseMult");
        dimage->photonCount.noise_mult = atof(argv[++i]);
    }
    if (!strncasecmp(argv[i], "-rhoVrhoG", 9))
    {
        checkArguments(1, i, argc, "-rhoVrhoG");
        dimage->photonCount.rhoVrhoG = atof(argv[++i]);
    }
    if (!strncasecmp(argv[i], "-photHDF", 8))
    {
        dimage->photonCount.writeHDF = 1;
#endif
    }