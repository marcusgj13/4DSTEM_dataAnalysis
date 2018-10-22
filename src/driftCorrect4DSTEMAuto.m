function [ s4DSTEM ] = driftCorrect4DSTEMAuto( s4DSTEM, offset )
% This function takes in a struct containig coordinates of single electron
% events from a 4DSTEM experiment and using centre of mass calculates the
% drift over the course of a scan.
%
% Marcus Gallagher-Jones 2018/10/20
% UCLA Department of Chemistry and Biochemistry
% marcusgj@chem.ucla.edu
% modification of original drift correct to now automatically identify the
% location of hot pixels and the central beam.

%% Extract dimensions of the 4D stack, identify the location of hot pixels
%  and the primary beam, and set up a place holder for individual frames of
%  diffraction.

xx = s4DSTEM.cubeSize(1);
yy = s4DSTEM.cubeSize(2);

tempIm = s4DSTEM.CBEDelectrons;
maxVal = max(tempIm(:));
badInds = find(tempIm >= maxVal/20);
tempIm(badInds) = 0;
smoothedIm = smooth2d(tempIm,0.3333);
[centY, centX] = find(smoothedIm == max(smoothedIm(:)));

s4DSTEM.shiftedElectrons = cell(s4DSTEM.cubeSize(1),s4DSTEM.cubeSize(2));
s4DSTEM.shiftedCBEDelectrons = zeros(s4DSTEM.cubeSize(3:4));
shift_array = zeros(2,xx);
%% Loop through the first dimension and find relative offsets of the images
%  based on the centre of mass of the central beam

for i = 1:xx
    img = makeClusterImage_strip(s4DSTEM,i,1:yy,2,0);
    [~,shifts] = COM_Align_subregion(img,centX - offset,centX + offset,centY - offset,centY + offset);
    shift_array(:,i) = shifts;
end

disp('=============shifts calculated for first dimension===========')
%% Shift all images and record the shifted indeces

for ii = 1:yy
    for jj = 1:xx
        img = makeClusterImage(s4DSTEM,jj,ii,0,0);
        img(badInds) = 0;
        img = circshift(img,[shift_array(1,jj),shift_array(2,jj)]);
        indsElectrons = find(img > 0);
        s4DSTEM.shiftedCBEDelectrons(indsElectrons) = ...
            s4DSTEM.shiftedCBEDelectrons(indsElectrons) + img(indsElectrons);
        [xe,ye] = ind2sub(s4DSTEM.cubeSize(3:4),indsElectrons);
        s4DSTEM.shiftedElectrons{jj,ii} = ...
            [xe ye img(indsElectrons) zeros(numel(xe),3)];
    end
    
end

shift_arrayFirst(1,:) = shift_array(1,:) + centY;
shift_arrayFirst(2,:) = shift_array(2,:) + centX;

disp('==============Shifting in first dimension complete===========')
%% repeat process but in second dimension
tempIm = s4DSTEM.shiftedCBEDelectrons;
smoothedIm = smooth2d(tempIm,0.3333);
[centY, centX] = find(smoothedIm == max(smoothedIm(:)));

s4DSTEM.electrons = s4DSTEM.shiftedElectrons;
s4DSTEM.shiftedCBEDelectrons = zeros(s4DSTEM.cubeSize(3:4));
shift_array = zeros(2,yy);

for i = 1:yy
    img = makeClusterImage_strip(s4DSTEM,1:xx,i,2,0);
    [~,shifts] = COM_Align_subregion(img,centX - offset,centX + offset,centY - offset,centY + offset);
    shift_array(:,i) = shifts;
end

disp('=============shifts calculated for second dimension===========')

for ii = 1:yy
    for jj = 1:xx
        img = makeClusterImage(s4DSTEM,jj,ii,0,0);
        img = circshift(img,[shift_array(1,ii),shift_array(2,ii)]);
        indsElectrons = find(img > 0);
        s4DSTEM.shiftedCBEDelectrons(indsElectrons) = ...
            s4DSTEM.shiftedCBEDelectrons(indsElectrons) + img(indsElectrons);
        [xe,ye] = ind2sub(s4DSTEM.cubeSize(3:4),indsElectrons);
        s4DSTEM.shiftedElectrons{jj,ii} = ...
            [xe ye img(indsElectrons) zeros(numel(xe),3)];
    end
    
end
disp('==============Shifting in second dimension complete===========')
shift_arraySecond(1,:) = shift_array(1,:) + centY;
shift_arraySecond(2,:) = shift_array(2,:) + centX;

%% While loop for minimizing the shifts

while mean(abs(shift_array(:))) > 2
    %% Loop through the first dimension and find relative offsets of the images
    %  based on the centre of mass of the central beam
    tempIm = s4DSTEM.shiftedCBEDelectrons;
    smoothedIm = smooth2d(tempIm,0.3333);
    [centY, centX] = find(smoothedIm == max(smoothedIm(:)));
    
    s4DSTEM.electrons = s4DSTEM.shiftedElectrons;
    s4DSTEM.shiftedCBEDelectrons = zeros(s4DSTEM.cubeSize(3:4));
    shift_array = zeros(2,xx);
    
    for i = 1:xx
        img = makeClusterImage_strip(s4DSTEM,i,1:yy,2,0);
        [~,shifts] = COM_Align_subregion(img,centX - offset,centX + offset,centY - offset,centY + offset);
        shift_array(:,i) = shifts;
    end
    
    disp('=============shifts calculated for first dimension===========')
    %% Shift all images and record the shifted indeces
    
    for ii = 1:yy
        for jj = 1:xx
            img = makeClusterImage(s4DSTEM,jj,ii,0,0);
            img = circshift(img,[shift_array(1,jj),shift_array(2,jj)]);
            indsElectrons = find(img > 0);
            s4DSTEM.shiftedCBEDelectrons(indsElectrons) = ...
                s4DSTEM.shiftedCBEDelectrons(indsElectrons) + img(indsElectrons);
            [xe,ye] = ind2sub(s4DSTEM.cubeSize(3:4),indsElectrons);
            s4DSTEM.shiftedElectrons{jj,ii} = ...
                [xe ye img(indsElectrons) zeros(numel(xe),3)];
        end
        
    end
    
    shift_arrayFirst(1,:) = shift_array(1,:) + centY;
    shift_arrayFirst(2,:) = shift_array(2,:) + centX;
    
    disp('==============Shifting in first dimension complete===========')
    %% repeat process but in second dimension
    tempIm = s4DSTEM.shiftedCBEDelectrons;
    smoothedIm = smooth2d(tempIm,0.3333);
    [centY, centX] = find(smoothedIm == max(smoothedIm(:)));
    
    s4DSTEM.electrons = s4DSTEM.shiftedElectrons;
    s4DSTEM.shiftedCBEDelectrons = zeros(s4DSTEM.cubeSize(3:4));
    shift_array = zeros(2,yy);
    
    for i = 1:yy
        img = makeClusterImage_strip(s4DSTEM,1:xx,i,2,0);
        [~,shifts] = COM_Align_subregion(img,centX - offset,centX + offset,centY - offset,centY + offset);
        shift_array(:,i) = shifts;
    end
    
    disp('=============shifts calculated for second dimension===========')
    
    for ii = 1:yy
        for jj = 1:xx
            img = makeClusterImage(s4DSTEM,jj,ii,0,0);
            img = circshift(img,[shift_array(1,ii),shift_array(2,ii)]);
            indsElectrons = find(img > 0);
            s4DSTEM.shiftedCBEDelectrons(indsElectrons) = ...
                s4DSTEM.shiftedCBEDelectrons(indsElectrons) + img(indsElectrons);
            [xe,ye] = ind2sub(s4DSTEM.cubeSize(3:4),indsElectrons);
            s4DSTEM.shiftedElectrons{jj,ii} = ...
                [xe ye img(indsElectrons) zeros(numel(xe),3)];
        end
        
    end
    disp('==============Shifting in second dimension complete===========')
    shift_arraySecond(1,:) = shift_array(1,:) + centY;
    shift_arraySecond(2,:) = shift_array(2,:) + centX;
end
disp('============== Shifting has stagnated ===========')
imagesc(sqrt(s4DSTEM.shiftedCBEDelectrons)); axis image; colormap gray;
hold
scatter(shift_arrayFirst(2,:),shift_arrayFirst(1,:),50,'r','filled');
scatter(shift_arraySecond(2,:),shift_arraySecond(1,:),50,'b','filled');
end

%% Subfunctions

function [shifted_image,shifts] = COM_Align_subregion(image_in,startX,endX,startY,endY)
%% Function for aligning an image stack based on the COM of a subregion
%  Within the image

[bigY, bigX, nIm] = size(image_in); % Dimensions of the input image
shifted_image = zeros(bigY,bigX,nIm); % Holder for aligned images

bigImCentY = ceil(bigY/2);
bigImCentX = ceil(bigX/2);
subCentY = ceil((startY+endY)/2);
subCentX = ceil((startX+endX)/2);

centOffsetY = bigImCentY - subCentY;
centOffsetX = bigImCentX - subCentX;

small_stack = image_in(startY:endY,startX:endX,:); % subregion of the image for alignment


[~,shifts] = COM_Align(small_stack); % align subregion first to get offsets
shifts(:,1) = shifts(:,1) + centOffsetY;
shifts(:,2) = shifts(:,2) + centOffsetX;

for ii = 1:nIm % Apply offsets to full image
    shifted_image(:,:,ii) = circshift(image_in(:,:,ii),[shifts(ii,1) shifts(ii,2)]);
end
end

function [shifted_array,shift] = COM_Align(imageN)
% Function for calculating alignment of a series of images based on
% The Images Centre of Mass

[xx1,yy1,zz1] = size(imageN); %Image dimensions

shifted_array = zeros(xx1,yy1,zz1); % Holder for the new images
centres = zeros(1,2,zz1); % Holder for the correct centres
shift = zeros(zz1,2); % holder for shifts
im_centrey = floor(xx1/2)+ 1;
im_centrex = floor(yy1/2) + 1;

for kk = 1:zz1
    A1 = smooth2d(abs(imageN(:,:,kk)),0.25);
    A2 = imageN(:,:,kk);
    C=cellfun(@(n) 1:n, num2cell(size(A1)),'uniformoutput',0);
    [C{:}]=ndgrid(C{:});
    C=cellfun(@(x) x(:), C,'uniformoutput',0);% create two vectors defining location of all pixels
    C=[C{:}];                                 % in x and y axis
    
    CenterOfMass = A1(:).'*C/sum(A1(:),'double'); % Calculate the image's Centre of mass
    
    centres(1,:,kk) = CenterOfMass;
    shift(kk,1) = round(im_centrex - centres(1,1,kk)); % Determine the amount to shift each image
    shift(kk,2) = round(im_centrey - centres(1,2,kk));
    shifted_array(:,:,kk) = circshift(A2,[shift(kk,1) shift(kk,2)]);
    
end
end

