function [s4DSTEM] = counting4DSTEM_03cluster( ...
    cube,CBEDbg,CBEDsub,coefsOrThresh,CBEDMean)
%% This function takes a Gatan 4DSTEM image stack with dimensions
%  [x_probe y_probe q_x q_y] and deconstructs images into electron
%  events using hybrid counting and stores this in a struct where each cell
%  is an Nx6 array with fields [xe ye I 0 0 0]
%
% Marcus Gallagher-Jones 
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu
%
% Colin Ophus
% National Center for Electron Microscopy, LBNL
% 2018/19/09


tic

numberSigmaThresh = 5;

% Take a Gatan 4D STEM cube with
% [x_probe y_probe q_x q_y] and find individual electron events
% place into a cell array of dimension
% [length(x_probe) length(y_probe)] where each cell is an
% N x 6 array, with fields
% [xe ye 0 0 0 0]


% Store everything in struct s4DSTEM.
if length(coefsOrThresh) > 1
s4DSTEM.coefs = coefsOrThresh;
s4DSTEM.CBEDbg = CBEDbg;
s4DSTEM.CBEDsub = CBEDsub;
s4DSTEM.CBEDMean = CBEDMean;
s4DSTEM.threshCluster = coefsOrThresh(2) * numberSigmaThresh + coefsOrThresh(3);
else
    s4DSTEM.threshCluster = coefsOrThresh;
end


% NN computation
s4DSTEM.cubeSize = size(cube);
indsNN = zeros(s4DSTEM.cubeSize(3)*s4DSTEM.cubeSize(4),8);
[ya,xa] = meshgrid(1:s4DSTEM.cubeSize(4),1:s4DSTEM.cubeSize(3));
xa = xa(:);
ya = ya(:);
indsNN(:,1) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-2,s4DSTEM.cubeSize(3))+1,mod(ya-2,s4DSTEM.cubeSize(4))+1);
indsNN(:,2) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-1,s4DSTEM.cubeSize(3))+1,mod(ya-2,s4DSTEM.cubeSize(4))+1);
indsNN(:,3) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-0,s4DSTEM.cubeSize(3))+1,mod(ya-2,s4DSTEM.cubeSize(4))+1);
indsNN(:,4) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-2,s4DSTEM.cubeSize(3))+1,mod(ya-1,s4DSTEM.cubeSize(4))+1);
indsNN(:,5) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-0,s4DSTEM.cubeSize(3))+1,mod(ya-1,s4DSTEM.cubeSize(4))+1);
indsNN(:,6) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-2,s4DSTEM.cubeSize(3))+1,mod(ya-0,s4DSTEM.cubeSize(4))+1);
indsNN(:,7) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-1,s4DSTEM.cubeSize(3))+1,mod(ya-0,s4DSTEM.cubeSize(4))+1);
indsNN(:,8) = sub2ind(s4DSTEM.cubeSize(3:4),...
    mod(xa-0,s4DSTEM.cubeSize(3))+1,mod(ya-0,s4DSTEM.cubeSize(4))+1);



% Perform hybrid electron counting
s4DSTEM.electrons = cell(s4DSTEM.cubeSize(1),s4DSTEM.cubeSize(2));
s4DSTEM.CBEDelectrons = zeros(s4DSTEM.cubeSize(3:4));
progressbar(0,2);
for ax = 1:s4DSTEM.cubeSize(1)
    for ay = 1:s4DSTEM.cubeSize(2)
        I = double(squeeze(cube(ax,ay,:,:))) - CBEDbg;
       
        indsCand = find(I > s4DSTEM.threshCluster);
        
        indsElectrons = indsCand(I(indsCand) ...
            > max(I(indsNN(indsCand,:)),[],2));
    
         I = floor((double(squeeze(cube(ax,ay,:,:))) - CBEDbg...
            )./(s4DSTEM.threshCluster));
       
        s4DSTEM.CBEDelectrons(indsElectrons) = ...
            s4DSTEM.CBEDelectrons(indsElectrons) + I(indsElectrons);

        [xe,ye] = ind2sub(s4DSTEM.cubeSize(3:4),indsElectrons);
        s4DSTEM.electrons{ax,ay} = ...
            [xe ye  I(indsElectrons) zeros(numel(xe),3)];
    end
    
    comp = ax / s4DSTEM.cubeSize(1);
    progressbar(comp,2);
end




% Testing plots
figure(1)
clf
imagesc(s4DSTEM.CBEDelectrons)
axis equal off
colormap(gray(256))
set(gca,'position',[0 0 1 1])

toc
end