function out = makeCircleMask(radius,imgSizex,imgSizey,centrex,centrey)

if nargin <= 3
    centrex = floor(imgSizex/2) + 1;
    centrey = floor(imgSizey/2) + 1;
end

[xx, yy] = meshgrid(1:imgSizex,1:1:imgSizey);
R = sqrt((xx-centrex).^2 + (yy-centrey).^2);
out = R<=radius;