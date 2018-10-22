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

    


function [shifted_array,shift] = COM_Align(imageN)
% Function for calculating alignment of a series of images based on 
% The Images Centre of Mass

[xx,yy,zz] = size(imageN); %Image dimensions

shifted_array = zeros(xx,yy,zz); % Holder for the new images
centres = zeros(1,2,zz); % Holder for the correct centres
shift = zeros(zz,2); % holder for shifts
im_centrey = floor(xx/2)+ 1;
im_centrex = floor(yy/2) + 1;

for kk = 1:zz;
    A1 = smooth2d(abs(imageN(:,:,kk)),0.25);
    A2 = imageN(:,:,kk);
    C=cellfun(@(n) 1:n, num2cell(size(A1)),'uniformoutput',0);
    [C{:}]=ndgrid(C{:});
    C=cellfun(@(x) x(:), C,'uniformoutput',0);% create two vectors defining location of all pixels
    C=[C{:}];                                 % in x and y axis
    
    CenterOfMass = A1(:).'*C/sum(A1(:),'double'); % Calculate the image's Centre of mass
%     display(CenterOfMass);
%     display(kk);
    centres(1,:,kk) = CenterOfMass;
    shift(kk,1) = round(im_centrex - centres(1,1,kk)); % Determine the amount to shift each image
    shift(kk,2) = round(im_centrey - centres(1,2,kk));
    shifted_array(:,:,kk) = circshift(A2,[shift(kk,1) shift(kk,2)]);

end
end
end
