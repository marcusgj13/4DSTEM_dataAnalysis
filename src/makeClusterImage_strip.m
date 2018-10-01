function [Icluster] = makeClusterImage_strip(s4DSTEM,xRange,yRange,sigma,shifted)
%% This function takes in a 4DSTEM counting struct and will out put an image
% using counts from a specified range of indeces within the struct.

% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu

Icluster = zeros(s4DSTEM.cubeSize(3),s4DSTEM.cubeSize(4));
for ii = xRange
    for jj = yRange
        if shifted == 1
            xyInds = s4DSTEM.shiftedElectrons{ii,jj};
        else
            xyInds = s4DSTEM.electrons{ii,jj};
        end
        inds = sub2ind(s4DSTEM.cubeSize(3:4),xyInds(:,1),xyInds(:,2));
        temp = zeros(s4DSTEM.cubeSize(3),s4DSTEM.cubeSize(4));
        temp(inds) = xyInds(:,3);
        Icluster = Icluster + temp;
    end
end


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