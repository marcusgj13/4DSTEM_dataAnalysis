function [s4DSTEM] = dm4toCounts4DSTEM(filename, xx, yy, zz1, zz2)
%% Function for converting a dm4 stack of 4DSTEM images and converting it 
%  to a struct containing lists of hybrid counts. This code will perform
%  backrgound estimation and subtraction, fitting of detector gaussian
%  background and counting threshold and finally conversion to hybrid
%  counts.
% WARNING - This code requires a machine with fairly high RAM to handle the
% large size of the dm4 stacks
%
% filename: .dm4 file to be loaded into memory 
% xx: scanX dimension
% yy: scanY dimension
% zz1: first image dimension (y in matlab)
% zz2: second image dimension (x in matlab)

% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu

%% Load and reshape the 4D stack of images

dm4Reader4D(filename,1); % creates DM4_tags.txt file to get scan parameters
outputData = dm4Reader4D(filename); % load data into memory
imageStack = outputData.cube; % output is a 1D vector need to reshape based on scan dimensions
imageStack = reshape(imageStack,[xx yy zz1 zz2]);

clear outputData

%% Calculate background image using median filter and subtract from all datasets
CBEDMean = squeeze(mean(mean(imageStack,1),2));
[CBEDsub,CBEDbg] = counting4DSTEM_01bg(CBEDMean);

%% Calculate threshold for considering electron events then perform counting

[coefs] = counting4DSTEM_02measThresh(imageStack,CBEDbg);
[s4DSTEM]= counting4DSTEM_03cluster(imageStack,CBEDbg,CBEDsub,coefs,CBEDMean);