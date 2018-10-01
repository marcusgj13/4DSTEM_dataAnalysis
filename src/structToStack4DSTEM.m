function [ rStack ] = structToStack4DSTEM( s4DSTEM, binFac, radius )
% This function converts a struct containg hybrid counts to a stack of
% binned images for use in clustering. This function will also perform 
% masking of the central beam to remove its influence on the clustering
% result.
%
% Marcus Gallagher-Jones 2018/09/19
% UCLA Department of Chemistry and Biochemistry
% marcusgj@chem.ucla.edu

%% Extract image dimensions and define size of new images and make holders
dimensions = s4DSTEM.cubeSize;

xx = dimensions(1);
yy = dimensions(2);
binY = round(dimensions(3)/binFac);
binX = round(dimensions(4)/binFac);
binStack = zeros(binY,binX,xx,yy);

%% Create a quick test image to roughly figure out the location of the central beam

testIm = imresize(s4DSTEM.shiftedCBEDelectrons,1/binFac);
[beamCentY, beamCentX] = find(smooth2d(testIm,0.5) == max(max(smooth2d(testIm,0.5))));

%% create a circular mask to remove the central beam from binned images
circMask = makeCircleMask(radius,binX,binY,beamCentX,beamCentY);

%% Loop through all of the images in the struct and resize them
disp('======== beginning to resize images ========')

totIm = xx*yy;
count = 1;
for ii = 1:xx
    for jj = 1:yy
        tempIm = imresize(makeClusterImage(s4DSTEM,ii,jj,0,1),1/binFac);
        tempIm(circMask == 1) = 0;
        binStack(:,:,ii,jj) = tempIm;
        count = count + 1;
        
        if mod(count,100) == 0
            disp([int2str(count) ' out of ' int2str(totIm)...
                ' patterns resized.'])
        end
    end
end

rStack = reshape(binStack,[binY binX xx*yy]);
clear binStack

end

