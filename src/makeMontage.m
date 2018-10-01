function montage_array = makeMontage(stackIn,numImX)
%% creates a single matrix containing a montage of images from an input
%  stack (stack_in). num_im_x defines the number of images along the first 
%  dimension
%
% Marcus Gallagher-Jones 2018/09/19
% UCLA Department of Chemistry and Biochemistry
% marcusgj@chem.ucla.edu

[xx,yy,zz] = size(stackIn); % get array dimensions

startX = 1; % initialise the loop for x and y
startY = 1;
xEnd = xx; % step size in x and y
yEnd = yy;

xLength = numImX*xx; % total length in the first dimension
yLength = numImX*yy;

for i = 1:zz
    % place image into matrix
    montage_array(startX:xEnd,startY:yEnd) = stackIn(:,:,i);
    startX = startX + xx; % move to next panel
    xEnd = xEnd + xx;
    if startX >= xLength % signal to move along the second dimension
        startX = 1;
        xEnd = xx;
        startY = startY + yy;
        yEnd = yEnd + yy;
    end
end