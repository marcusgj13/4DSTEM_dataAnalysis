# 4DSTEM_dataAnalysis
A set of matlab functions for preprocessing and analyzing data from 4DSTEM experiments.
Author: Marcus Gallagher-Jones  
Institution: UCLA Department of chemistry and biochemistry  
Email: marcusgj13@gmail.com 

## Introduction
This set of functions were used to analyze low-dose 4DSTEM diffraction data from nanocrystals and determine the orientation of nm sized regions of the crystals. The associated publication can be found here (add hyperlink at some point). These functions can be used to first extract meaningful data from raw diffraction patterns and then perform unsupervised clustering to group regions of simmilar diffraction. At some point I may convert these scripts over to python also.

### Dependencies
Running these functions requires Matlab 2013A or later and the image processing and curve fitting toolboxes.

### Instructions
#### Initial preprocessing
The first step in analyzing 4DSTEM data is to reduce the raw images (in .dm4 format) into lists of hybrid counts, _KxKy_ location on the detector, and their associated intensity. This is performed by the __dm4toCounts4DSTEM.m__ function which calls the subfunctions __counting4DSTEM_01bg.m__, __counting4DSTEM_02measThresh.m__ and __counting4DSTEM_03cluster.m__. 
The code can be called as follows
```matlab
s4DSTEM = dm4toCounts4DSTEM('path/to/file.dm4', scanX, scanY, zz1, zz2);
```
Where scanX and ScanY are the number of scan points and zz1 and zz2 are the two principle dimensions of the diffraction image in pixels (1792 and 1920 for a K2-IS detector). The output is then a struct (s4DSTEM) containing the following categories:

| Struct entry | Description |
| --- | --- |
| coeffs | The coefficients of the gaussian fit to dark noise used to create counting threshold. |
| CBEDbg | Differential darkcurrent offset of the detector calculated by median filter. |
| CBEDsub |  Mean diffraction image calculated from the 4DSTEM image stack after background subtraction. |
| CBEDMean| Mean diffraction image calculated from the 4DSTEM image stack before background subtraction. |
| threshCluster | Threshold used for determining hybrid counts. |
| coeffs | The coefficients of the gaussian fit to dark noise used to create counting threshold. |
| cubeSize | Dimensions of the original 4DSTEM data stack. |
| electrons | Cell containing lists of _KxKy_ count locations and intensity values for all images within the 4DSTEM stack. |
| CBEDelectrons | Diffraction pattern created by summing all true counts across all images. |

#### Shift correction
Following conversion to hybrid counts individual patterns will have to be corrected for x and y shift so that they all have a common centre. This is dependent on the realspace size of the scan. Larger scan dimensions (several microns) will have more severe beam shift and the correction may need to be run multiple times. This is handled by the function __driftCorrect4DSTEM.m__ which aligns all images to a single point using the COM of the central disk. After running the struct will be updated to contain two new entries, __shiftedElectrons__ and __shiftedCBEDElectrons__, which are the new _KxKy_ location of the counts and the summed pattern respectively:
```matlab
s4DSTEM = driftCorrect4DSTEM(s4DSTEM, badInds, centreX, centreY, offset, startShifted);
```
__badInds__ is a vector containing the _XY_ location of hot pixels within the image, __centreX and Y__ are the pixel coordinates of the centre of the image, __offset__ is the size in pixels of the region around the centre to crop for COM alignment of the central disk, __startShifted__ indicates wether or not to use the shifted electron positions as a start point if running the script multiple times. All of the data provided has already gone through these two steps.

#### Preparing data for clustering
Clustering is performed on a reduced for of the images. This is to both improve SNR and also to reduce computation times. This is performed by the function __structToStack4DSTEM.m__ which will bin all images within the struct by a user defined ammount and also mask out the central beam. The rationale for this is that the central beam contains the most intensity and as such will heavily bias the outcome of clustering, rather than the presence/absence of Bragg peaks. The function is run as follows:
```matlab
rStack = structToStack4DSTEM(s4DSTEM, binFactor, radius);
```
__binFactor__ is the ammount you wish to reduce the image, for example to reduce the image to 1/8th its current size binFactor should be equal to 8. __radius__ is the radius of a circular mask used to remove the central disk from all binned images. Setting this too large may have strange effects on the clustering.

#### Unsupervised clustering of 4DSTEM data
Clustering is performed via K-means with initial clusters being assigned by the Kmeans++ algorithm and simmilarity determined via Euclidean distance. Optionally _K_ can also be determined by G-means. Note that G-means has been modified such that 80% of the clustered data must be found to be gaussian by the Anderson-Darling statistic instead of 100% to account for clusters over free space being close to unity due to the low background noise. This algorithm is not determenistic and as such will give slightly different results with different runs if initialized with different cluster centres. The code is run as follows:
```matlab
[ clusteredInds, mindisFunc, wcss, meanImages, minShifts, initialCentres ] = ...
    KMeansPP4DSTEM( stackIn, numSets, numIterations,doMean, doShifts, initialize, scanX, scanY);
```
__numIterations__ is the total number of iterations to run, generally the algorithm converges before 20 - 30 iterations have run so I set this to 100. __doMean__ determines whether or not to initialise the algorithm with an image that is the mean of several images nearby the intial image or just start with a single image. Generally performance is best when using a local mean as the information content of a single image is fairly low. __doShifts__ determines wether or not to perform image alignment during the clustering, either by cross-correlation or circshifting. This is generally not necessary for diffraction data. __initialize__ should be a list of centre coordinates, this is useful when trying to run the code reproducably from the same starting point.

Following clustering it is recommended to run __regularizeClusters4DSTEM.m__ which will remove any unphysical single pixel clusters resulting from misassignment during the K-means clustering.

#### Orientation assignment by library matching
The first step in orientation assignment is to generate a library of NBED patterns from a known structure. This is performed by the __makeTilts4DSTEM.m__ function. First the coordinates of all atoms in a small nanocrystal should be generated. The function __makeProteinCell_QyN9.m__ will generate atomic coordinates and lattice parameters for a crystal structure of [this](http://www.rcsb.org/structure/6AXZ) peptide structre. The function can be modified for other structures by providing a vector containing atom positions in fractional coordinates. The library is generated as follows:
```matlab
[tiltLibrary, tilts] = makeTilts4DSTEM(tiltRange, doGPU, outputName);
```
__tiltRange__ is a vector of angles in degrees, __doGPU__ determines whether or not to run the __STEM_tilts__ function on the GPU or not, __outputName__ is the name to give to the .mat file that contains simulated images and their tilt angles. Note: if you can run on GPU, run on GPU otherwise the library generation is __VERY__ slow.

Once the library has been simulated it can be passed, along with the mean diffraction patterns determined by the K-means clustering , to the __orientationSearch4DSTEM.m__ function. This function will match the diffraction pattern to a simulated diffraction pattern by comparing the RMSD of normalised intensity values at all potential Bragg peak locations. It is run as follows:
```matlab
[bestMatch, bestErr, simStruct, expStruct] = orientationSearch4DSTEM(...
    simTiltIms, tilts, meanImages, centres);
```
__centres__ represents the location of all peaks, and potential peaks within the diffraction pattern. This can be calculated from a binned version of shiftedCBEDelectrons using __RealspaceLattice01.m__ or any other function you might have for peak finding (template matching also works well). __bestMatch__ is a vector containing the index of the simulated pattern that as the lowest RMSD (stored in __bestErr__) to the experimental pattern. __bestAngles__ is the _XY_ tilt relative to some fixed orientation of the best matched diffraction pattern, i.e. the orientation of the lattice represented by the mean diffraction pattern that was fit. This can then be used to replot the cluster map in terms of orientation using the function __plotTiltsRGB.m__.
