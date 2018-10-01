function [tiltLibrary, tilts] = makeTilts4DSTEM(tiltRange, doGPU, outputName)
%% This function will generate a library of simulated NBED patterns to use
%  with the library mathching algorithm in orientationSearch4DSTEM.m to
%  determine the orientation of the experimental clusters. This can either
%  be run on CPU or GPU (recommended).
%
%  tiltRange: range of tilts in the x direction (the range for y should be
%  set in STEMtilt(GPU)_QYN9.m. This is mostly for memory limitation
%  purposes.
%  doGPU: whether or not to perform simulations with GPU
%  outputName: name of the .mat file to save the library of patterns and
%  the tilts

% Marcus Gallagher-Jones
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu


counter = 1;

for i = tiltRange
    [atoms,cellDim] = makeProteinCell_QyN9();
    
    if doGPU == 1
        [emdSTEM] = STEMtilts_gpu_QYN9(atoms, cellDim, i, 1);
    else
        [emdSTEM] = STEMtilts_QYN9(atoms, cellDim, i, 1);
    end
    
    for jj = 1
        for ii = 1:41
            SimPattern = rot90(fftshift(squeeze(squeeze(emdSTEM.data(:,:,40,1,jj,ii)))));
            tiltLibrary(:,:,:,counter) = SimPattern;
            tilts(1,counter) = emdSTEM.xTiltArray(jj);
            tilts(2,counter) = emdSTEM.yTiltArray(ii);
            counter = counter + 1;
        end
    end
    disp([num2str(counter) ' out of 1681 patterns simulated'])
    close all
end

save(outputName,'tiltLibrary', 'tilts','-v7.3');
