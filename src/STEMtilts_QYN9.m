function [emdSTEM] = STEMtilts_QYN9(atoms,cellDim,numFP)
%% This function simulates NBED diffraction patterns over a range
%  of thickness that would be collected from a 4DSTEM experiment. It
%  requires a known structure which is provided by the two inputs, atoms
%  and cellDim, and also knowledge of the microscope parameters. Data is
%  returned as a struct and the individual patterns are found under
%  emdSTEM.data.
%
% Marcus Gallagher-Jones
% Department of Chemistry and Biochemistry, UCLA
% marcusgj@chem.ucla.edu
%
% Colin Ophus
% National Center for Electron Microscopy, LBNL
% 2018/19/09

tic

numUC = [2 2];  % Tiling of UCs (to prevent probe boundary artifacts)
thicknessOutput = 100:100:4000;
f_plotPot = 1;   % Plot total projected potential (for testing)
% also shows simulated probe positions!

defocusArray = 0;  % vector of probe defocus values in Angstroms
xTiltArray = xTiltArray*pi/180;
yTiltArray = (-4:0.2:4)*pi/180;

% debye waller coefs (in units of RMS atomic displacement in Angstroms)
u = ones(118,1)*0.1;

% simulation parameters
% Probe positions
Np = [1 1]*3;
xp = linspace(0,cellDim(1),Np(1)+1); xp(end) = [];
yp = linspace(0,cellDim(2),Np(2)+1); yp(end) = [];


potBound = 0.8;
alphaMax = 0.2/1000;  % illumination semiangle in rads

E0 = 300*10^3;  % Microscope voltage in volts
pSize = 0.2;
subSlice = 128 * 3;  % full cell
flagKeepPot = 0;


if nargin == 2
    numFP = 1;
end

% Generate simulation supercell
[yUC,xUC] = meshgrid(0:(numUC(2)-1),0:(numUC(1)-1));
tile = [xUC(:) yUC(:) zeros(numel(xUC),2)];
[indTile,indAtoms] = meshgrid(1:size(tile,1),1:size(atoms,1));
atomsAll = atoms(indAtoms,:) + tile(indTile,:);

% Scale atomic positions by cellDim
atomsAll(:,1) = atomsAll(:,1)*cellDim(1);
atomsAll(:,2) = atomsAll(:,2)*cellDim(2);
atomsAll(:,3) = atomsAll(:,3)*cellDim(3);

% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
lambda = h/sqrt(2*m*e*E0)/sqrt(1 + e*E0/2/m/c^2) * 10^10; % wavelength in A
s = (2*pi/lambda/E0)*(m*c^2+e*E0)/(2*m*c^2+e*E0);

% Make realspace coordinate systems
Nx = ceil(numUC(1)*cellDim(1)/pSize/4)*4;  % Make sure factor 4 pixels per cell
Ny = ceil(numUC(2)*cellDim(2)/pSize/4)*4;
xSize = numUC(1)*cellDim(1) / Nx;
ySize = numUC(2)*cellDim(2) / Ny;
xySize = [xSize ySize];


% find z plane value for sub slices
zPlanes = mod(round(atomsAll(:,3)/cellDim(3)*subSlice) - 1, subSlice) + 1;
zPlanesAll = unique(zPlanes);
dz = cellDim(3) / subSlice;
thicknessArray = max(round(thicknessOutput / dz),1);
thickness = thicknessArray * dz;

emdSTEM = struct();
emdSTEM.xySize = xySize;
emdSTEM.E0 = E0;
emdSTEM.lambda = lambda;
emdSTEM.sigma = s;
emdSTEM.numFP = numFP;
emdSTEM.dz = dz;
emdSTEM.xp = xp;
emdSTEM.yp = yp;
emdSTEM.thickness = thickness;
emdSTEM.numUC = numUC;
emdSTEM.defocusArray = defocusArray;
emdSTEM.xTiltArray = xTiltArray;
emdSTEM.yTiltArray = yTiltArray;
emdSTEM.atoms = atoms;
emdSTEM.cellDim = cellDim;

% Make Fourier coordinates
Lx = Nx*xSize;
Ly = Ny*ySize;
qx = circshift(((-Nx/2):(Nx/2-1))/Lx,[0 -Nx/2]);
qy = circshift(((-Ny/2):(Ny/2-1))/Ly,[0 -Ny/2]);
[qya, qxa] = meshgrid(qy,qx);

q2 = qxa.*qxa + qya.*qya;
q1 = sqrt(q2);

xOutput = [(1:(Nx/4)) ((1-Nx/4):0)+Nx];
yOutput = [(1:(Ny/4)) ((1-Ny/4):0)+Ny];
emdSTEM.qxOutput = qx(xOutput);
emdSTEM.qyOutput = qy(yOutput);

% Make propagators and anti aliasing aperture AA
dq = qx(2) - qx(1);
Adist = (max(qx)/2 - q1)/dq+.5;
AA = Adist;
AA(Adist>1) = 1;
AA(Adist<0) = 0;

% Propagator Array
Ntilt = [length(xTiltArray) length(yTiltArray)];
propArray = zeros(Nx,Ny,Ntilt(1),Ntilt(2));
for a0 = 1:Ntilt(1)
    xTilt = xTiltArray(a0);
    for a1 = 1:Ntilt(2)
        yTilt = yTiltArray(a1);
        propArray(:,:,a0,a1) = exp(-1i*pi*lambda*dz*q2 ...
            + 2*pi*1i*dz*qxa*xTilt ...
            + 2*pi*1i*dz*qya*yTilt).*AA;
    end
end

% Probe
qMax =  alphaMax / lambda; 
Adist = (qMax - q1)/dq+.5;
A = Adist;
A(Adist>1) = 1;
A(Adist<0) = 0;

% Construct projected potentials
xyLeng = ceil(potBound./xySize);
xvec = -xyLeng(1):xyLeng(1);
yvec = -xyLeng(2):xyLeng(2);
xr = xvec*xySize(1);
yr = yvec*xySize(2);

% Lookup table for atom types
atomTypes = unique(atoms(:,4));
potLookup = zeros(length(xvec),length(yvec),length(atomTypes));
uLookup = zeros(length(atomTypes),1);
for a0 = 1:length(atomTypes)
    potLookup(:,:,a0) = projPot(atomTypes(a0),xr,yr);
    if u(atomTypes(a0)) > 0
        uLookup(a0) = u(atomTypes(a0));
    else
        disp(['Warning, no RMS displacement given for atomic number ' ...
            num2str(atomTypes(a0)) '!  Setting u = 0.1 A.'])
        uLookup(a0) = 0.1;
    end
end


% Generate potentials
potAll = zeros(Nx,Ny,subSlice);
for a1 = 1:length(zPlanesAll)
    % subset of atoms on this slice
    aSub = atomsAll(zPlanesAll(a1) == zPlanes,:);
  
    % Generate slice potential
    for a2 = 1:size(aSub,1)
        [~,ind] = min(abs(atomTypes-aSub(a2,4)));
        x = mod(xvec+round((aSub(a2,1)+randn*uLookup(ind))/xySize(1)),Nx)+1;
        y = mod(yvec+round((aSub(a2,2)+randn*uLookup(ind))/xySize(2)),Ny)+1;
        potAll(x,y,a1) = potAll(x,y,a1) + potLookup(:,:,ind);
    end
end
if flagKeepPot == true
    emdSTEM.potAll = potAll;
end
trans = exp(1i*s*potAll);
if f_plotPot == 1
    potSum = sum(potAll,3);
end


% Main loops

emdSTEM.data = zeros( ...
    Nx/2,Ny/2,...
    length(thickness),...
    length(defocusArray),Ntilt(1),Ntilt(2));
emdSTEM.intData = zeros(thicknessArray(end),2);
psi = zeros(Nx,Ny,length(xp),length(yp),length(defocusArray),Ntilt(1),Ntilt(2));
for a0 = 1:numFP
    % Initialize all probes
    for a3 = 1:length(defocusArray)
        chiProbe = pi*lambda*q2*defocusArray(a3);
        
        for a1 = 1:length(xp)
            for a2 = 1:length(yp)
                probefft = exp(-1i*chiProbe ...
                    - 2*pi*1i*(qxa*(xp(a1)) ...
                    + qya*(yp(a2)))).*A;
                probefft = probefft ...
                    / sqrt(sum(sum(abs(probefft(:)).^2)));
                for a4 = 1:Ntilt(1)
                    for a5 = 1:Ntilt(2)
                        psi(:,:,a1,a2,a3,a4,a5) = probefft;
                    end
                end
            end
        end
    end
    
    % Propagate all probes through sample
    for a1 = 1:thicknessArray(end)

        indSlice = mod(a1-1, subSlice) + 1;
        
        % apply slice phase shift
        for a2 = 1:length(xp)
            for a3 = 1:length(yp)
                for a4 = 1:length(defocusArray)
                    for a5 = 1:Ntilt(1)
                        for a6 = 1:Ntilt(2)
                            psi(:,:,a2,a3,a4,a5,a6) = fft2(ifft2(...
                                psi(:,:,a2,a3,a4,a5,a6)) ...
                                .*trans(:,:,indSlice)) ...
                                .*propArray(:,:,a5,a6);
                        end
                    end
                end
            end
        end
        emdSTEM.intData(a1,:) = emdSTEM.intData(a1,:) ...
            + [dz*a1 ...
            sum(abs(psi(:).^2))/length(xp)/length(yp)/Ntilt(1)/Ntilt(2)];
        
        % Output results if needed
        [val,indOut] = min(abs(thicknessArray - a1));
        if val == 0            
            for a2 = 1:length(xp)
                for a3 = 1:length(yp)
                    for a4 = 1:length(defocusArray)
                        for a5 = 1:Ntilt(1)
                            for a6 = 1:Ntilt(2)
                                emdSTEM.data(:,:,indOut,a4,a5,a6) = ...
                                    emdSTEM.data(:,:,indOut,a4,a5,a6) ...
                                    + abs(psi(xOutput,yOutput,...
                                    a2,a3,a4,a5,a6)).^2;
                            end
                        end
                    end
                end
            end
        end
        
        % timing
        comp = (a1 / (thicknessArray(end)) ...
            + a0 - 1) / numFP;
        progressbar(comp,2);
    end
    
end
emdSTEM.data = emdSTEM.data / numFP / length(xp) / length(yp);

if f_plotPot == 1
    emdSTEM.potSum = potSum/numFP;
    
    figure(1)
    clf
    imagesc(potSum/numFP)
    hold on
    [yy,xx] = meshgrid(yp/xySize(2)+1,xp/xySize(1)+1);
    scatter(yy(:),xx(:),'r.')
    hold off
    axis equal off
    colormap(gray(256))
    set(gca,'position',[0 0 1 1])
end




toc
end
