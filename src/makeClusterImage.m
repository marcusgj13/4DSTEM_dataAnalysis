function [Icluster] = makeClusterImage(s4DSTEM,xInd,yInd,sigma,shifted)
%% This function takes in a 4DSTEM counting struct and will out put an image
% using counts from a specified index within the struct.

% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu


Icluster = zeros(s4DSTEM.cubeSize(3),s4DSTEM.cubeSize(4));
if shifted == 1
    xyInds = s4DSTEM.shiftedElectrons{xInd,yInd};
else
    xyInds = s4DSTEM.electrons{xInd,yInd};
end
inds = sub2ind(s4DSTEM.cubeSize(3:4),xyInds(:,1),xyInds(:,2));
Icluster(inds) = xyInds(:,3);

if nargin == 3
    sigma = 0;
end


if sigma > 0
    x = 1:s4DSTEM.cubeSize(3);
    y = 1:s4DSTEM.cubeSize(4);
    x = mod(x + s4DSTEM.cubeSize(3)/2,s4DSTEM.cubeSize(3)) - s4DSTEM.cubeSize(3)/2;
    y = mod(y + s4DSTEM.cubeSize(4)/2,s4DSTEM.cubeSize(4)) - s4DSTEM.cubeSize(4)/2;
    [ya,xa] = meshgrid(y,x);
    kernel = exp(-(xa.^2 + ya.^2) / (2*sigma^2));
    Icluster = ifft2(fft2(Icluster) .* fft2(kernel) ,'symmetric');
    Icluster = max(Icluster,0);
end

end