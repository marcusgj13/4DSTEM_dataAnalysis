function [CBEDsub,CBEDbg] = counting4DSTEM_01bg(CBEDmean)
%% This function takes in the mean diffraction pattern calculated from a 
%  4DSTEM image stack and calculates the detector dark current offset using
%  a median filter.
% 
% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu
%
% Colin Ophus
% National Center for Electron Microscopy, LBNL
% 2018/19/09
% Range to estimate background from along y direction, outside ADF aper

yRange = [];
xRange = [1200 1792];

imageSize= size(CBEDmean);
if ~isempty(yRange)
    CBEDbg = repmat(median( ...
        CBEDmean(:,yRange(1):yRange(2)),2),[1 imageSize(2)]);
else
    CBEDbg = repmat(median( ...
        CBEDmean(xRange(1):xRange(2),:),1),[imageSize(1) 1]);
end
% Subtract background from CBED image

CBEDsub = CBEDmean - CBEDbg;

quickImagePlot([CBEDmean-median(CBEDmean(:)) CBEDsub],[-0.5 0.5+2]);

end