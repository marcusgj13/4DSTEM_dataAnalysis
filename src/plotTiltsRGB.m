function [ imageRGB ] = plotTiltsRGB(indexPoints,angleIndices,indexErr)
%% Takes in a set of XY tilts observed in a 4DSTEM experiment and converts
%  them to RGB coordinates after normalization.
%  
% Marcus Gallagher-Jones 2018/09/19
% UCLA Department of Chemistry and Biochemistry
% marcusgj@chem.ucla.edu

% Normalize coordinates such that zero tilt is equal to a value of 0.5
scaleFac = 1;
bestX = angleIndices(1,:) - mean(angleIndices(1,:));
bestY = angleIndices(2,:) - mean(angleIndices(2,:));
maxValX = max(abs(bestX));
maxValY = max(abs(bestY));
bestX = (bestX - maxValX.*-1)./(2*maxValX);
bestY = (bestY - maxValY.*-1)./(2*maxValY);

% initialize empty vectors
xyOutput = [80 82];
Ir = zeros(xyOutput);
Ig = zeros(xyOutput);
Ib = zeros(xyOutput);

% mask from empty pixels
mask = ~(indexPoints == 1);
mask = reshape(mask, xyOutput);

numClasses = size(angleIndices,2);
for a0 = 2:numClasses
    Ir(indexPoints==a0) = bestY(a0);
    Ib(indexPoints==a0) = bestX(a0);
    Ig(indexPoints==a0) = indexErr(a0);
    
end

imageRGB = zeros(xyOutput(1),xyOutput(2),3);
imageRGB(:,:,1) = Ir.*mask.*scaleFac;
imageRGB(:,:,2) = Ig.*mask;
imageRGB(:,:,3) = Ib.*mask.*scaleFac;

figure; imagesc(imageRGB)
axis equal off
end

