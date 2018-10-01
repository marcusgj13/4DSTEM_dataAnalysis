function [ clusteredInds, mindisFunc, wcss, meanImages, minShifts, initialCentres ] = ...
    KMeansPP4DSTEM( stackIn, numSets, numIterations,doMean, doShifts, initialize, scanX, scanY)
% This function performs K means clustering on a stack of images to seperate
% them into K sets of simmilar images.
%   stackIn = stack of images to be sorted.
%   numSets = number of sets to split the data set into for clustering
%   numIterations = number of iterations of clustering to perform
%   doMean = initialize from the mean of a local cluster
%   doShifts = can be either 0,1 or 2. O is no shift, 1 is norm
%   Xcorrelation and 2 is by circshifting
%   initialize = provide a set of indices for preinitializing the algorithm
%
%
% Marcus Gallagher-Jones 2018/09/19
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu
%
% Modified 05/12/2017: include K++ initialization
% Modified 14/02/2018: include clustering by cross-correlation
% Modified 26/02/2018: include flag for gpu usage
% Modified 27/02/2018: include circshift into euclidean distance
% calculation
% Modified 28/02/2018 include preinitialization and saving initial cluster
% centres
%% reset the random number generator
rng('shuffle', 'twister');
tic
%% Read in the data and initialize starting parameters

if nargin < 3
    numIterations = 30;
end
num_rows = floor(sqrt(numSets));
g_count = 0;%gpuDeviceCount; % check for gpu

[xx, yy, zz] = size(stackIn); % Determine number of images
% holder to track which cluster a particular image is in
clusteredInds = zeros(1,zz);
% holder to track mean distance of images to cluster centres
if g_count >= 1
    disFunc = gpuArray(zeros(numSets,zz,'single'));
    mindisFunc = gpuArray(zeros(1,zz,'single'));
    meanImages = gpuArray(zeros(xx,yy,numSets,'single'));
    dist = gpuArray(zeros(1,zz,'single'));
    initialCentres = gpuArray(zeros(1,numSets,'single'));
    minShifts = gpuArray(zeros(2,zz,'single'));
else
    disFunc = zeros(numSets,zz,'single');
    mindisFunc = zeros(1,zz,'single');
    meanImages = zeros(xx,yy,numSets,'single');
    dist = zeros(1,zz,'single');
    initialCentres = zeros(1,numSets,'single');
    minShifts = zeros(2,zz,'single');
end

% within-cluster sum of squares - start with large value
wcss = 10e9;
% holder for mean images
iteration = 1;


%% Assign initial centres using the kmeans ++ algorithm described by Arthur
%  Vassilvitskii (2006)
if initialize == 0
    randomizedIndeces = randperm(zz);
    if doShifts == 1
        meanImages(:,:,1:numSets) = stackIn(:,:,randomizedIndeces(1:numSets));
    else
        meanImages(:,:,1:numSets) = stackIn(:,:,randomizedIndeces(1:numSets));
    end
    
    fprintf('================ Beginning assigning centres ==============\n')
    for cc = 2:numSets+1
        summation = 0;
        for ii = 1:zz
            [~, min_dist] = findNearestCluster(stackIn(:,:,ii),...
                meanImages(:,:,1:cc-1), 0);
            dist(ii) = min_dist;
            summation = summation + dist(ii);
        end
        
        summation = summation.*rand;
        for jj = 1:size(dist,2)
            summation = summation - dist(jj);
            if summation > 0
                continue
            end
            if doMean == 1
                if jj > size(dist,2) - 3
                    meanImages(:,:,cc-1) = nanmean(stackIn(:,:,jj-6:jj),3);
                elseif jj < 4
                    meanImages(:,:,cc-1) = nanmean(stackIn(:,:,jj:jj+6),3);
                else
                    meanImages(:,:,cc-1) = nanmean(stackIn(:,:,jj-3:jj+3),3);
                end
            else
                meanImages(:,:,cc-1) = stackIn(:,:,jj);
                initialCentres(cc-1) = jj;
            end
            break
        end
        
        
        
    end
else
    meanImages(:,:,1:numSets) = stackIn(:,:,initialize);
    initialCentres = initialize;
end
fprintf('================ Done assigning centres ==============\n')

%% Calculate starting distance functions and wcss

for jj = 1:zz
    [min_ind, min_dist, shifts] = findNearestCluster(stackIn(:,:,jj),...
        meanImages,doShifts);
    clusteredInds(jj) = min_ind;
    mindisFunc(jj) = min_dist;
    minShifts(:,jj) = shifts;
end

for ii = 1:numSets
    if size(clusteredInds(clusteredInds == ii),2) > 0
        clusterIms = stackIn(:,:,clusteredInds == ii);
        clusterShifts = minShifts(:,clusteredInds == ii);
        for jj = 1:size(clusterIms,3)
            clusterIms(:,:,jj) = circshift(clusterIms(:,:,jj),[clusterShifts(1,jj) clusterShifts(2,jj)]);
        end
        meanImages(:,:,ii) = squeeze(nanmean(clusterIms,3));
    end
end

old_wcss = wcss;
wcss = sum(mindisFunc);

%% Iterate until stagnation or desired number of iterations

while iteration < numIterations && old_wcss > wcss
    for jj = 1:zz
        [min_ind, min_dist, shifts] = findNearestCluster(stackIn(:,:,jj),...
            meanImages,doShifts);
        clusteredInds(jj) = min_ind;
        mindisFunc(jj) = min_dist;
        minShifts(:,jj) = shifts;
    end
    
    for ii = 1:numSets
        if size(clusteredInds(clusteredInds == ii),2) > 0
            clusterIms = stackIn(:,:,clusteredInds == ii);
            clusterShifts = minShifts(:,clusteredInds == ii);
            for jj = 1:size(clusterIms,3)
                clusterIms(:,:,jj) = circshift(clusterIms(:,:,jj),...
                    [clusterShifts(1,jj) clusterShifts(2,jj)]);
            end
            meanImages(:,:,ii) = squeeze(nanmean(clusterIms,3));
        end
    end
    
    % Just show the middle portion of the mean images
    subplot(1,2,1); imagesc(makeMontage(meanImages,num_rows)); 
    axis image off; colormap jet2black
    title('Cluster centres')
    subplot(1,2,2); imagesc(reshape(clusteredInds,[scanX scanY])); 
    axis image off; colormap jet2black
    title(['Iteration = ' int2str(iteration) ' clusters = ' int2str(numSets)])
    drawnow
  
    old_wcss = wcss;
    wcss = sum(mindisFunc);
    
    fprintf('============ Iteration %d complete ============\n', iteration);
    iteration = iteration + 1;
    
end
if g_count >= 1
    [wcss, minShifts, meanImages, mindisFunc, initialCentres] = ...
        gather(wcss, minShifts, meanImages, mindisFunc, initialCentres);
end
t = toc;
fprintf('=========== Total time elapsed = %f mins =========\n',t/60)
end
%% Subfunctions
%% Calculate the square distance between two points
function di = sqr_distance(p, cc)
di = sum((cc(:)-p(:)).^2);
end


%% Calculate relative offsets using circshift and euclidean distance

function [di, best_shifts] = circ_sqr_dist(p,cc,displacement)
best_shifts = [0 0];
di = 10000;

for ii = -displacement:displacement
    for jj = -displacement:displacement
        p_shift = circshift(p,[ii jj]);
        dist = sqr_distance(p_shift,cc);
        if dist < di
            di = dist;
            best_shifts(1) = ii;
            best_shifts(2) = jj;
        end
    end
    
end
end

%% Calculate a correlation coefficient and relative offsets for each image
%  relative to a cluster

function [di, shifts] = get_XCorr(p,cc, apMask)
p = padarray(p,[3 3]).*apMask;
cc = padarray(cc,[3 3]).*apMask;
xc = normxcorr2_general(cc,p);
di = max(xc(:));
[xCen, yCen] = find(xc == di,1,'first');
shifts(1) = xCen - size(cc,1);
shifts(2) = yCen - size(cc,2);
di = 1 - di; % This term is added so the overall function is still a minimization.
end

%% Find the nearest cluster centre
function [min_ind, min_dist, best_shifts] = findNearestCluster(p,...
    cluster_centres,do_shifts)

min_ind = 0;
best_shifts = zeros(2,1);
min_dist = 10000;
if do_shifts == 1
    apMask= abs(apMask);
end

if do_shifts == 2
    displacement = 6;
end

for dd = 1:size(cluster_centres,3)
    if do_shifts == 1
        [d, shifts] = get_XCorr(p,cluster_centres(:,:,dd),apMask);
        if min_dist > d
            min_dist = d;
            min_ind = dd;
            best_shifts = shifts;
        end
    elseif do_shifts == 2
        [d, shifts] = circ_sqr_dist(p,cluster_centres(:,:,dd), displacement);
        if min_dist > d
            min_dist = d;
            min_ind = dd;
            best_shifts = shifts;
        end
    else
        d = sqr_distance(p,cluster_centres(:,:,dd));
        if min_dist > d
            min_dist = d;
            min_ind = dd;
        end
    end
end
end