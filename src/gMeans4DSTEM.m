function [ clusteredInds, mindisFunc, wcss, meanImages, shifts ] = gMeans4DSTEM( stackIn, doXcorr, scanX, scanY )
%% This function takes in a set of 4DSTEM images (stackIn)and performs K-means
%  clustering using a combination of Kmeans ++ to properly initialise and G
%  means clustering to 'learn' the correct K from the data assuming that
%  data is clustered in a gaussian way about the cluster centres.
% 
% Marcus Gallagher-Jones 2018/09/19
% Department of Chemistry and Biochemistry, UCLA

numSets = 1; % number of sets of K to initialize with
numIterations = 50; % number of iterations of Kmeans to run
exit_status = 0;
while exit_status < 1
    prevNumsets = numSets;
    [ clusteredInds, mindisFunc, wcss, meanImages, shifts ] = ...
        KMeansPP4DSTEM( stackIn, numSets, numIterations,1,doXcorr, scanX, scanY);
    pTot = zeros(1,prevNumsets);
    for i = 1:numSets
        
        try
            if size(mindisFunc(clusteredInds == i),2) >= 10
                set = mindisFunc(clusteredInds == i);
                set = set./std(set);
                set = set - mean(set);
                [h, p] = adtest(set,'Alpha', 0.001,'MCTol',0.001);
                pTot(i) = p;
                display(p)
                if h == 1
                    numSets = numSets + 1;
                end
            else
                pTot(i) = 1;
            end
            
            
        catch
            display('Too few datapoints in a cluster exiting.')
            exit_status = 1;
        end
    end
    percentGauss = size(find(pTot >= 0.001),2)/prevNumsets;
    display(percentGauss);
    if percentGauss >= 0.8
        exit_status = 1;
    end
    if numSets == prevNumsets
        break
    end
end

[ clusteredInds, mindisFunc, wcss, meanImages, shifts ] = ...
    KMeansPP4DSTEM( stackIn, prevNumsets, numIterations,1,doXcorr, scanX, scanY);

end

