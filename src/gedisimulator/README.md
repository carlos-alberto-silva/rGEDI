# README #

### What is this repository for? ###

This is a set of programs to simulate large-footprint full-waveform lidar from airborne lidar and to process it and perform various other tasks. It is described and validated in:

[Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner, J.R. and Dubayah, R., 2019. The GEDI simulator: A large‐footprint waveform lidar simulator for calibration and validation of spaceborne missions. Earth and Space Science.](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018EA000506)


A google group has been set up for asking questions. It can be found [here](https://groups.google.com/g/gedisim).


The programs are:


**gediRat**: simulates GEDI waveforms from ALS .las files and outputs ASCII or HDF5 waveforms.

**gediMetric**: processes full-waveform data (LVIS or simulated GEDI) and outputs metrics.

**mapLidar**: produces geotiffs from .las files of different properties. Also produces the bounding boxes of .las files for use in overlapLasFiles.csh.

**lasPoints**: outputs .pts files from .las files for selected areas.

**collocateWaves**: collocates GEDI or LVIS to ALS data, following [Blair and Hofton (1999)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/1999GL010484).

**addNoiseHDF**: Reads waveform data from HDF5 files and adds noise of a chosen level.


To find options, type the above command with "-help". These are explained in full detail later. 


There are some C and bash shells to control the above:

**gediRatList.csh**: batch processes gediRat.

**filtForR.csh**: converts the output .txt files into .csv files for reading in to R.

**overlapLasFiles.cs**h: determines which .las files are needed for a given simulation.

**orbitTracks.bash**: produces lists of footprints from GEDI orbital simulations.

**listALS.csh**: produces the lists of .las files needed to read multiple files.

The other .c files are either small test programs in the development of the above or are libraries called by the above. Important libraries are:

**gediIO.c**: contains the GEDI simulator functions.

### How do I get set up? ###

There are three ways to install the code.

* Singularity container: Simplest but needs root access to set up
* Compile from source: Do this if you are confident compiling C code
* Compilation script: Automatically does the above. Do this if you don't have root access and are not confident compiling C code

These three methods are listed below. Note that the compilation script will create a directory structure to install the code in to and modify your .bashrc file to point to these, so it is not recommended if you have a particular way you like your system set up.


##### Singularity container
The simplest way to get set up is to use the script provided to build all programs within a [**singularity container**](https://www.sylabs.io/docs/). Download this [script](https://bitbucket.org/StevenHancock/gedisimulator/src/master/makeSingularity.txt), then build (in Linux) by the command below. **Note** that this requires root access. Singularity cannot be run on a Mac natively and must be built and run within a [virtual machine](https://www.sylabs.io/guides/2.6/user-guide/installation.html#install-on-mac).

    sudo singularity build gediSingularity makeSingularity.txt

This will compile all programs within a singularity, called *gediSingularity* (the name can be changed to whatever the user desires). All the above programs can be called, for example:

    singularity exec --bind $dataDir gediSingularity gediMetric -help

Where *dataDir* is a disk location that you want the singularity to have access to. Each disk area the singularity needs access to must be bound in this way.



#### Compile from source

Alternatively to the singularity container, the programs can be compiled from source. Clone the repository and point to its location with the environment variable:

    GEDIRAT_ROOT


All programs depend on these libraries:

* Gnu Scientific Library
* Geotiff
* HDF5
* GDAL
* [Minpack (Levenberg-Maruqardt)](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html)


Point to these with the environment variables:

    GSL_ROOT
    HDF5_LIB
    CMPFIT_ROOT

Once they're installed, it also requires two libraries without package managers:

* [C-tools](https://bitbucket.org/StevenHancock/tools)
* [libClidar](https://bitbucket.org/StevenHancock/libclidar)

Clone these from bitbucket and point to their locations with the environment variables:

    HANCOCKTOOLS_ROOT
    LIBCLIDAR_ROOT

To compile, type:

  **make THIS=gediRat**

Make sure that ~/bin/$ARCH exists and is in the path, where $ARCH is the result of `uname -m`.

  **make THIS=gediRat install**

Replace "**gediRat**" with each of the commands above to compile and install.


Make sure that all **.csh** and **.bash** files are also in your path.


### Compilation script

There is a bash script that will create a directory structure, clone the necessary libraries that do not have package managers and compile the script. *Note* that the package managed libraries are still required and these should be installed first (Gnu Scientific Library, Geotiff, HDF5, GDAL).

Download and run the script by typing the following commands in your terminal (Unix or Linux):

    wget https://bitbucket.org/StevenHancock/gedisimulator/src/master/installGedi.bash
    chmod +x installGedi.bash
    installGedi.bash

The will create the directory:

   ***$HOME/src***

To clone the source code to and will compile the executables to:

   ***$HOME/bin***

It will add that directory to your PATH so that your system will be able to find the commands. Afterwards you can test with:

   collocateWaves -help

Which, if it has worked, will print out the options for the collocateWaves tool.



# Tool  operation #

## gediRat ##

Program to create GEDI waveforms from ALS las or pts files. laz not yet supported. Data is output either as ASCII files or as a HDF5 file, both of which can be ready by gediMetric below.



#### Usage example

gediRat reads from ALS data in .las format and outputs waveforms in either ASCII or HDF5 format. To read data from a single las file (``file.las''), simulate a single footprint at coordinates **lon, lat**, and write the results to an ASCII file, use the following command:

    gediRat -input file.las -coord $lon $lat -output waveform.txt

To read data from multiple las files, perform a grid of simulations and write the output to a HDF5 file, use the following command:

    gediRat -inList alsList.txt -output waveforms.h5 -hdf -ground -step 25 -gridBound $minX $maxX $minY $maxY

Where ``alsList.txt'' is an ASCII file containing a list of absolute filenames (including path) of las files. This will make a grid of simulations between minX,minY and maxX,maxY with 25 m steps. The output will be HDF5 format and the ground portion of the waveforms will be labelled.




##### Input output filenames and format
    -input name;     lasfile input filename
    -inList list;    input file list (ASCII file) for multiple files
    -output name;    output filename
    -ground;         record separate ground and canopy waveforms
    -hdf;            write output as HDF5. Best with gridded or list of coords
    -ascii;          write output as ASCII (default). Good for quick tests
    -waveID id;      supply a waveID to pass to the output (only for single footprints)

##### Single footprint, list of footprints, or grid of footprints
    -coord lon lat;  footprint coordinate in same system as lasfile
    -listCoord name; list of coordinates
    -gridBound minX maxX minY maxY;     make a grid of waveforms in this box
    -gridStep res;   grid step size

##### Lidar characteristics. Defaults are expected GEDI values.
    -pSigma sig;     set Gaussian pulse width as 1 sigma
    -pFWHM fhwm;     set Gaussian pulse width as FWHM in ns. This is the outgoing laser pulse, which will be twice the detected pulse width
    -readPulse file; read pulse shape and width from a file instead of making Gaussian
    -fSigma sig;     set footprint width
    -wavefront file; read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma
    -res res;        range resolution of waveform digitisation to output, in units of ALS data
    -LVIS;           use LVIS pulse length, sigma=6.25m
    -topHat;         use a top hat wavefront
    -sideLobe;       use side lobes
    -lobeAng ang;    lobe axis azimuth


##### Input data quality filters
    -checkCover;     check that the footprint is covered by ALS data. Do not output if not
    -maxScanAng ang; maximum scan angle, degrees
    -decimate x;     probability of accepting an ALS beam


##### Computational speed options
    -pBuff s;        point reading buffer size in Gbytes
    -maxBins;        Optional: for HDF5, limit number of bins to save trimming.
    -countOnly;      only use count method
    -pulseAfter;     apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)
    -pulseBefore;    apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed
    -noNorm;         don't normalise for ALS density

##### Octree
    -noOctree;       do not use an octree
    -octLevels n;    number of octree levels to use
    -nOctPix n;      number of octree pixels along a side for the top level


##### Using full-waveform input data (not tested)
    -decon;          deconvolve
    -indDecon;       deconvolve individual beams
    -readWave;       read full-waveform where available


##### Miscellaneous
    -listFiles;      list files. Do not read them
    -keepOld;        do not overwrite old files, if they exist
    -useShadow;      account for shadowing in discrete return data through voxelisation
    -polyGround;     find mean ground elevation and slope through fitting a polynomial
    -nnGround;       find mean ground elevation and slope through nearest neighbour


### A note on ALS data

The **-ground** flag in gediRat uses the ALS ground classification to separate the contributions to the waveform into ground and canopy. This weighting between ground and canopy is then used in gediMetric to estimate the ALS derived canopy cover. Therefore, the tolerance used when classifying the ground points in the ALS data may affect the canopy cover estimate returned by gediMetric if too tight a tolerance is used. The ground classification tolerance should be enough to account for any topographic variations within a footprint as well as any ALS noise (particularly at the edge of flightlines). For operational use, gediRat is run on ALS data that has been classified with a 60 cm tolerance.



## gediMetric ##

Program to process large-footprint lidar data (real or simulated) and produce standard waveform metrics. It can add noise to simulations and alter pulse shapes (increase length only). It reads either ASCII or HDF5 files created by gediRat, or can read LVIS data in either HDF5 or .lgw format. It will be updated to read GEDI data when that is available. Take care when reading ASCII data as some options are mutually exclusive (different gediRat options can change the column order). This outputs an ASCII file with the first row defining the contents of each column. Output variable names are defined below.


#### Usage example

gediMetric reads waveform files in either ASCII, HDF5 (GEDI L1B, LVIS or gediRat output) or lgw (LVIS binary) format. It needs to be told what format the input is in. An example command on a HDF5 GEDI file (either simulated by gediRat or L1B) is:

    gediMetric -input waveforms.h5 -readHDFgedi -ground -varScale 3.5 -sWidth 0.8 -rhRes 2 -laiRes 5

That will use a noise threshold of the mean plus 3.5\*standard deviation, a smoothing width of 0.8 m and then output RH metrics in 2% intervals and the LAI profile (PAVD) in 5 m intervals.


##### Input output
    -input name;      waveform  input filename
    -outRoot name;    output filename root
    -inList list;     input file list for multiple files
    -writeFit;        write fitted waveform
    -writeGauss;      write Gaussian parameters
    -readBinLVIS;     input is an LVIS binary file
    -readHDFlvis;     read LVIS HDF5 input
    -readHDFgedi;     read GEDI simulator HDF5 input
    -level2 name;     level2 filename for LVIS ZG
    -bounds minX minY maxX maxY;    only analyse data within bounds

##### Switches
    -ground;          read true ground from file
    -useInt;          use discrete intensity instead of count
    -useFrac;         use fractional hits rather than counts
    -rhRes r;         percentage energy resolution of RH metrics
    -laiRes res;      LAI profile resolution in metres. Default 10 m.
    -laiH h;          height to calculate LAI to
    -noRHgauss;       do not fit Gaussians
    -gTol tol;        ALS ground tolerance. Used to calculate slope.
    -fhdHistRes res;  waveform intensity resolution to use when calculating FHD from histograms
    -forcePsigma;     do not read pulse sigma from file
    -bayesGround;     use Bayseian ground finding
    -dontTrustGround; don't trust ground in waveforms, if included
    -noRoundCoord;    do not round up coords when outputting

##### Adding noise:
    -dcBias n;        mean noise level
    -nSig sig;        noise sigma
    -seed n;          random number seed
    -hNoise n;        hard threshold noise as a fraction of integral
    -linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover
    -linkFsig sig;    footprint width to use when calculating and applying signal noise
    -linkPsig sig;    pulse width to use when calculating and applying signal noise
    -trueSig sig;     true sigma of background noise
    -nSkew s;         skewness of  background noise
    -nPeriodAmp a;    periodic noise amplitude
    -nPeriodOm p;     periodic noise wavelength
    -bitRate n;       digitisation bit rate
    -maxDN max;       maximum DN
    -renoise;         remove noise from truth before applying new noise level
    -newPsig sig;     new value for pulse width, when lengthening pulse
    -oldPsig sig;     old value for pulse width if not defined in waveform file, when lengthening pulse
    -addDrift xi;     apply detector background drift
    -missGround;      assume ground is missed to assess RH metrics
    -minGap gap;      delete signal beneath min detectable gap fraction

##### Photon counting
    -photonCount;     output point cloud from photon counting
    -nPhotons n;      mean number of photons
    -photonWind x;    window length for photon counting search, metres
    -noiseMult x;     noise multiplier for photon-counting

##### Denoising:
    -meanN n;         mean noise level, if using a predefined mean level
    -thresh n;        noise threshold, if using a predefined noise threshold
    -varNoise;        use a variable noise threshold
    -varScale x;      variable noise threshold scale (multiple of stdev above mean to set threshold)
    -statsLen len;    length to calculate noise stats over for varNoise
    -noiseTrack;      use noise tracking
    -sWidth sig;      smoothing width, after densoising
    -psWidth sigma;   smoothing width, before denoising
    -msWidth sig;     smoothing width, after noise stats, before denoising
    -preMatchF;       matched filter before denoising
    -postMatchF;      matched filter after denoising
    -pFile file;      read pulse file, for deconvolution and matched filters
    -gWidth sig;      Gaussian parameter selection smoothing width
    -minGsig sig;     minimum Gaussian sigma to fit
    -minWidth n;      minimum feature width in bins
    -medNoise;        use median stats rather than mean
    -varDrift;        correct detector drift with variable factor
    -driftFac xi;     fix drift with constant drift factor
    -rhoG rho;        ground reflectance
    -rhoC rho;        canopy reflectance
    -pSigma sig;      pulse width to smooth by if using Gaussian pulse
    -gold;            deconvolve with Gold's method
    -deconTol;        deconvolution tolerance

### gediMetric variable names ###
Note that some metrics are "true" and will not be available to GEDI. They are included to assess errors and sensitivities.


##### Metrics available to GEDI
    gHeight - ground elevation (m) from Gaussian fitting
    gSlope - ground slope (degrees) from Gaussian fitting
    maxGround - ground elevation (m) from lowest maximum
    inflGround - ground elevation (m) from inflection points.
    signal top - elevation of first point above noise (may include noise tracking).
    signal bottom - elevation of last return above noise (may include noise tracking).
    cover - canopy cover (fraction) from area of Gaussian fitted ground. Uses rho_v=0.57 and rho_g=0.4.
    leading edge ext - leading edge extent (m), from Lefksy et al (2007).
    trailing edge extent - trailing edge extent (m), from Lefksy et al (2007).
    rhGauss 0-100 - RH metrics, 0%-100%, using ground from Gaussian fitting (m).
    rhMax 0-100 - RH metrics, 0%-100%, using ground from lowest maximum (m).
    rhInfl 0-100 - RH metrics, 0%-100%, using ground from inflection points (m).
    gaussHalfCov - canopy cover (fraction) from double the energy beneath the Gaussian ground. Uses rho_v=0.57 and rho_g=0.4.
    maxHalfCov - canopy cover (fraction) from double the energy beneath the lowest maximum ground. Uses rho_v=0.57 and rho_g=0.4.
    infHalfCov - canopy cover (fraction) from double the energy beneath the inflection point ground. Uses rho_v=0.57 and rho_g=0.4.
    bayHalfCov - canopy cover (fraction) from double the energy beneath the experimental "Bayesian" ground. Uses rho_v=0.57 and rho_g=0.4.
    lon - footprint centre longitude in projection of ALS data (m).
    lat - footprint centre latitude in projection of ALS data (m).
    waveEnergy - total energy within waveform (will be 1 scaled by noise for simulations).
    blairSense - Blair's sensitivity metric. Canopy cover at which this SNR would have90% chance of detecting ground (does not account for rho_v/rho_g).
    FHD - Foliage height diversity
    niM2 - Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2 (unpublished)
    niM2.1 - Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2.1 (unpublished)
    gLAI0t10 - LAI at each height band (0-10m in this case) using the ground Gaussian to isolate canopy returns
    hgLAI0t10 - LAI at each height band (0-10m in this case) using the energy beneath the lowest Gaussian, reflected, to isolate canopy returns
    hiLAI0t10 - LAI at each height band (0-10m in this case) using the energy beneath the lowest inflection point, reflected, to isolate canopy returns
    hmLAI0t10 - LAI at each height band (0-10m in this case) using the energy beneath the lowest maximum, reflected, to isolate canopy returns

##### Metrics unavailable to GEDI
    wave ID - waveform label, relates to plot name and footprint number.
    true ground - ground elevation (m) from ALS. Centre of gravity of ground points within footprint
    true top - elevation of highest point of waveform (m), RH99.9, without noise. Includes pulse blurring.
    ground slope - effective ground slope (degrees), from width of ground return. Includes roughness.
    ALS cover - canopy cover (fraction) from ALS data. Uses rho_v=0.57 and rho_g=0.4.
    rhReal 0-100 - RH metrics, 0%-100%, using "true" ground from ALS data (m).
    groundOverlap - fraction of ground return overlapping with canopy return. A measure of understorey.
    groundMin - depth of minimum between ground and canopy return. A measure of understorey.
    groundInfl - d2y/dx2 of inflection point between ground and canopy return. A measure of understorey.
    pointDense - average ALS point density within GEDI footprint.
    beamDense - average ALS beam density within GEDI footprint.
    tLAI0t10 - LAI at each height band (0-10m in this case) using the ALS ground estimate to isolate canopy returns

##### System settings
    pSigma - GEDI system pulse width, sigma (m).
    fSigma - GEDI footprint width, sigma (m).
    linkM - link margin if noise is added (db).
    linkCov - canopy cover at which the above link margin is true (fraction).
    filename - name of input waveform file.




##### Signal processing description


##### Gaussian fitting
Used for "gHeight", "rhGauss" and "gaussHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75) and Gaussians fitted with Levenberg-Marquardt optimisation. The centre of the lowest Gaussian containing at least 0.5% of the waveform energy is selected as the ground.

##### Maximum
Used for "maxGround", "rhMax" and "maxHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75). The lowest maximum is taken as the ground.

##### Inflection points
Used for "inflGround", "rhInfl" and "inflHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75). The centre of gravity between the lowest two inflection points is taken as the ground.

##### Half covers
Used for "gaussHalfCov", "maxHalfCov" and "inflHalfCov". Sum energy beneath estimated ground position. Double that is the ground energy. Calculate canopy cover, correcting for rho_v and rho_g.

cover=Ecan/(Ecan+Eg*rho_v/rho_g)

Where Ecan is the canopy energy, Eg is the ground energy, rho_v is the vegetation reflectance and rho_g is the ground reflectance.


##### Edge extents
These are described in:

Lefsky, Michael A., Michael Keller, Yong Pang, Plinio B. De Camargo, and Maria O. Hunter. "Revised method for forest canopy height estimation from Geoscience Laser Altimeter System waveforms." Journal of Applied Remote Sensing 1, no. 1 (2007): 013537-013537.




## collocateWaves ##

Uses the correlation method in Blair and Hofton (1999) to colocate a large-footprint lidar dataset with a small-footprint, discrete-return dataset. Note that it requires the full-waveform LVIS or GEDI data, which is contained in the L1B files. It uses the Pearson correlation to find the best affine transformation (x and y only, or x, y and z) and footprint size needed to align a large-footprint dataset with a small-footprint dataset. It has three potential modes of operation.

* It can test a grid of affine transformations and give the correlation for every point (for a single footprint width), as used in Blair and Hofton (1999)
* It can use a simplex to move along the error surface and find the optimum transformation and footprint size. Note that initial location needs to be within around 20 m of the true location for their to be a sufficient gradient on the error surface.
* It can use a hybrid of the two, testing every step on a coarse grid, then setting a simplex off from the location of maximum correlation. This allows a rapid assessment to get within 20 m of the true value, then uses the simplex to find the offsets precisely.

It requires the large-footprint system pulse shape (either a Gaussian width or a file containing range and intensity) and the EPSG codes for the two datasets. It reads ALS data in .las format and can read either simulated HDF5 files from gediRat or LVIS in HDF5 or lgw format. A reader for GEDI data will be added once that data is available. 

If the full grid is used, it outputs an ASCII file with the correlation for each x, y and z offset. If it uses a simplex it outputs the single optimum offset.


#### Usage example

The following example is currently the most efficient for finding the offset between GEDI and ALS data and outputting simulations of GEDI from the ALS aligned with GEDI.

    collocateWaves -listALS alsList.txt -gedi waveforms.h5 -readHDFgedi -aEPSG 32622 -solveCofG -geoError 30 5 -fixFsig -writeWaves simulated.h5 -minDense 3 -minSense 0.9


Where ``alsList.txt'' is a list of ALS filenames with absolute path, wavefroms.h5 is the GEDI L1B or simulated file, ``32622'' is the EPSG code of the ALS data and ``simulated.h5'' is the simulated waveforms output filename. This will only use ALS data with at least 3 beams per square metre and only use GEDI waveforms with at least 90% beam sensitivity.


#### Input-output options
    -output name;     output filename
    -listAls list;    input file list for multiple als files
    -als file;        input als file
    -gedi file;       single input GEDI/LVIS L1B file
    -listGedi file;   list of multiple GEDI/LVIS files
    -readHDFgedi;     read GEDI HDF5 input (default)
    -lgw;             LVIS is in lgw (default is GEDI hdf5)
    -readHDFlvis;     read GEDI HDF5 input (default is GEDI hdf5)
    -bounds minX minY maxX maxY;    bounds to use, specified in ALS projection
    -leaveEmpty;      exit if there are no usable footprints
    -lEPSG epsg;      LVIS projection
    -aEPSG epsg;      ALS projection

#### Grid mode operation
    -maxShift x;      grid mode, horizontal distance to search over
    -step x;          grid mode, horizontal step size
    -maxVshift x;     grid or geoError mode, vertical distance to search over
    -vStep z;         grid or geoError mode, vertical step size

#### Optimiser mode operation
    -simplex;         use simplex optimisation rather than doing the full bullseye plot
    -anneal;          use simulated annealing optimisation
    -fixFsig;         fix fSigma in simplex
    -geoError expError correlDist;   rapid geolocation, using expected geolocation error and correlation distance. Vertical shifts must be separatley defined
    -quickGeo;        perform rapid geolocation using default error values
    -optTol x;        tolerance for optimisation
    -maxIter n;       maximum number of iterations
    -writeSimProg;    write progress of simplex to output
    -writeWaves name; write out final waveforms as HDF5 when using simplex
    -nTriesAnneal n;  how many points do we try before stepping?
    -itersFixedT n;   how many iterations for each T?
    -kAnneal x;       Boltzmann constant for annealing
    -tInitial x;A     initial annealing temperature
    -muAnneal x;      damping factor for temperature
    -tMinAnneal x;    minimum annealing temperature

#### Initial estimates. Will search around this point
    -hOffset dx dy;   centre of horizontal offsets
    -offset z;        vertical datum offset

#### Waveform characteristics
    -fSigma x;        footprint width, sigma in metres
    -pSigma x;        Gaussian pulse length, sigma in metres
    -readPulse file;  pulse shape, if not Gaussian

#### Filters for input data
    -minSense x;      minimum LVIS/GEDI beam sensitivity to accept
    -maxZen zen;      maximum LVIS/GEDI zenith angle to use, degrees
    -maxScanAng ang;  maximum ALS scan angle, degrees
    -minDense x;      minimum ALS beam density to accept
    -decimate f;      decimate ALS point cloud by a factor, to save RAM
    -noFilt;          don't filter outliers from correlation (default)
    -filtOut s;       filter outliers from correlation stats along with an optional sigma
    -smooth sig;      smooth both waves before comparing
    -checkCov;        check ALS coverage and remove any with footprints with less than 2/3 coverage
    -median;          use the median correlation when optimising, rather than the mean (default)

#### GEDI beam selection
    -beamList 11111111; 0/1 for whether or not to use beams 1-8 on GEDI
    -skipBeams n;     list of GEDI beam numbers to skip. No spaces between (eg 123)
    -readBeams n;     list of GEDI beam numbers to read. No spaces between (eg 123)

#### Simulator settings. For simulator validation only
    -noNorm;          don't correct sims for ALS densiy variations
    -norm;            correct sims for ALS densiy variations
    -allSimMeth;      use all simulation methods
    -pulseBefore;     apply pulse shape before binning to prevent aliasing

#### Octree to speed searching of ALS data. Not fully operational
    -noOctree;       do not use an octree
    -octLevels n;    number of octree levels to use
    -nOctPix n;      number of octree pixels along a side for the top level


## plotWaveComparison.py ##

Plots waveforms to compare ALS derived simulations to real GEDI data. It is meant to take the output from collocateWaves' -writeWaves option and compare to the original GEDI data. If the collocation has been successful, the two will match.

## mapLidar ##
Generates a geotiff from las file properties, combining multiple files. Can also print a list of file bounds or calculate beam and point density.


#### Options
    -input name;     lasfile input filename
    -output name;    output filename
    -inList list;    input file list for multiple files
    -res res;        image resolution, in metres
    -bounds minX minY maxX maxY;     user defined image bounds
    -float;          output as float
    -height;         draw height image
    -cover;          draw canopy cover map
    -noInt;          no image
    -findDens;       find point and footprint density
    -epsg n;         geolocation code if not read from file
    -writeBound n;   write file bounds to a file
    -pBuff s;        point reading buffer size in Gbytes
    -printNpoint;    print number of points in each file


## lasPoints ##
Extracts a point cloud as a pts for a bounding box within a collection of las files.


### Contribution guidelines ###

Please talk to svenhancock@gmail.com to suggest edits.


### License ###

Gnu Public License

### Who do I talk to? ###

Questions can be posted on this [Google group](https://groups.google.com/g/gedisim) and answers found.

svenhancock@gmail.com
