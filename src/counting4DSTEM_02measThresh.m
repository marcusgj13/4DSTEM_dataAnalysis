function [coefs] = counting4DSTEM_02measThresh(cube,CBEDbg)
%% This function fits a gaussian background to the intensity distribution
%  of all pixels within the 4DSTEM image stack to determine threshold
%  levels for hybrid counting.
%
% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu
%
% Colin Ophus
% National Center for Electron Microscopy, LBNL
% 2018/19/09

skip = 10; % skip this many entries in cube
ints = double(cube(1:skip:end,1:skip:end,:,:));
for a0 = 1:size(ints,1)
    for a1 = 1:size(ints,2)
        ints(a0,a1,:,:) = squeeze(ints(a0,a1,:,:)) - CBEDbg;
    end
end
ints = ints(:);


numberBins = [256 32];
includeThresh = 0.999;

% coefs
fun = @(c,x) c(1)*exp((-1/2/c(2)^2) ...
    *((x-c(3)).^2));
options = optimset('TolX',1e-5,'TolFun',1e-5,...
    'MaxFunEvals',10^4);

% histogram range
range = [min(ints) max(ints)];
x = linspace(range(1),range(2),numberBins(1));
dx = x(2) - x(1);
counts = histc(ints,x-dx/2);

% filter input data
sig = cumsum(counts);
sig = sig / sig(end);
[~,inds] = min(abs(sig - includeThresh));
intensityThresh = x(inds);

% Repeat histogram measurement
ints = ints(ints<intensityThresh);

% histogram range
range = [min(ints) max(ints)];
x = linspace(range(1),range(2),numberBins(2))';
dx = x(2) - x(1);
counts = histc(ints,x-dx/2);

% Fit histogram
sigma0 = sqrt(sum(x.^2 .* counts) / sum(counts));
coefs = [max(counts) sigma0 0];
[~,coefs] = evalc(['lsqcurvefit(fun,coefs,'  ...
    'x,counts,[],[],options)']);

% plot results
figure(1)
clf
scatter(x,counts,'k.','sizedata',80)
hold on
countsFit = fun(coefs,x);
plot(x,countsFit,'linewidth',2,'color','r')
hold off



end