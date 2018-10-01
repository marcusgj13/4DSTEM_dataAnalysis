function newClusterIm = regularizeClusters4DSTEM(clusterIm,kernalSize)
% This function takes in a cluster map from Kmeans and regularizes it using
% a square kernal to get rid of spuriously misassigned data points, e.g.
% patterns that have been assigned to the wrong cluster causing a big jump
% in orientation.
%
% Marcus Gallagher-Jones 2018/09/19
% UCLA Department of Chemistry and Biochemistry
% marcusgj@chem.ucla.edu
%% Extract the coordinates of the image and create a place holder for the new image
[xx, yy] = size(clusterIm);
newClusterIm = zeros(xx,yy);
offset = floor(kernalSize/2);
centrePixel = round((kernalSize^2)/2);
majority = floor((kernalSize.^2)/2);

%% Perform the regularization
for i = 2:xx-1
    for j = 2:yy-1
        currPixel = clusterIm(i,j);
        neighbours = clusterIm(i-offset:i+offset,j-offset:j+offset);
        neighbours = neighbours(:);
        neighbours(centrePixel) = [];
        mde = mode(neighbours(:));
        if currPixel ~= mde
            if size(find(neighbours == mde),1) > majority
                currPixel = mde;
            end
        end
        newClusterIm(i,j) = currPixel;
    end
end
newClusterIm(newClusterIm == 0) = clusterIm(newClusterIm == 0);