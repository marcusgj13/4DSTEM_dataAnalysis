function [bestMatch, bestErr, simStruct, expStruct] = orientationSearch4DSTEM(...
    simTiltIms, meanImages, centres)
% Function for loading a library of simulated 4DSTEM patterns at a range of
% tilts that have been scaled to match experimental patterns derived from
% clustering of 4DSTEM data.
%
% simTiltIms: A four dimensional stack of simulated NBED patterns from
% makeTilts4DSTEM.m. The first two dimensions are the image size the third
% dimension a range of thicknesses and the fourth the range of tilts.
% meanImages: Experimental patterns that come from KMeansPP4DSTEM.m or
% gMeans4DSTEM.m
% centres: Location of Bragg peaks within the mean CBED pattern, excluding
% the central beam.
% 
% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu

%% The library should first be filtered to reduce noise that is an
%  artifact of the simulation and then cropped.
numSim = size(simTiltIms,4);

% Convert centres to integers
centres = round(centres);

% Set up mask for each potential bragg spot
for i = 1:size(centres,1)
    peakMask(:,:,i) = makeCircleMask(8,240,224,centres(i,2),centres(i,1));
end
%% Set up structs to store the positions of all the peaks and their
% integrated intensities (need to do some scaling here too)

simStruct = struct();
expStruct = struct();

thickness = 100:100:6000;

%% Identify the positions of all the peaks and store them in struct
fprintf('====== Identifying peaks in simulated patterns =======\n')
for ii = 1:numSim
    for jj = 1:60 % For all thicknesses
    temp = squeeze(simTiltIms(:,:,jj,ii));
    
    % place holders for position and intensity to go into struct
    intensity = zeros(1,size(centres,2));
    % Loop over all identified bragg peaks and calculate intensity
    for kk = 1:size(centres,1)
        intensity(kk) = sum(temp(peakMask(:,:,kk) == 1));   
    end
    
    % add everything to struct
    simStruct.thickness{ii,jj,:} = thickness(jj);
    simStruct.tiltX{ii,jj} = tilts(1,ii);
    simStruct.tiltY{ii,jj} = tilts(2,ii);
    simStruct.I{ii,jj,:} = intensity;
    end
end

fprintf('====== Done identifying peaks =======\n')
%% Identify bragg peaks in experimental patterns in the same way
fprintf('====== Identifying peaks in experimental patterns =======\n')
for ii = 1:size(meanImages,3)
    temp = meanImages(:,:,ii);
    % Top hat filter reduces the diffuse intensity close to the primary
    % beam
    temp = imtophat(temp,strel('disk',10));
    
    % place holders for position and intensity to go into struct
    intensity = zeros(1,size(centres,2));
    % Loop over all identified bragg peaks and calculate intensity
    for kk = 1:size(centres,1)
        intensity(kk) = sum(temp(peakMask(:,:,kk) == 1));   
    end
    
  
    expStruct.I{ii,:} = intensity;
end

fprintf('====== Done identifying peaks =======\n')
%% Perform allignment by iteratively calculating the distance between 
%  Intensity of all peaks in simulated and experimental patterns.
%  Needs terms for background and scaling for intensity
bestErr = zeros(1,size(meanImages,3));
bestMatch = zeros(1,size(meanImages,3));

allErr = zeros(size(meanImages,3),numSim,length(thickness));

fprintf('====== Starting Peak matching algorithm =======\n')
for ii = 1:size(meanImages,3)
    expI = expStruct.I{ii,:};
    % Make intensities the ratio of the strongest spot
    expI = expI./max(expI(:));
    
    for jj = 1:num_sim
        for kk = 1:length(thickness)
            sim_I = simStruct.I{jj,kk,:};
            sim_I = sim_I./max(sim_I(:));
            allErr(ii,jj,kk) = sum(sqrt((expI - sim_I).^2))/length(expI);
        end
    end
end
fprintf('====== Finding the best matches ======\n')

for ii = 1:size(meanImages,3)
    bestErr(ii) = min(min(allErr(ii,:,:)));
    [xx, yy] = find(squeeze(allErr(ii,:,:)) == bestErr(ii));
    bestMatch(1,ii) = xx;
    bestMatch(2,ii) = yy;
   
end
fprintf('====== Done! =======\n')
    