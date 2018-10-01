function [I] = quickImagePlot(I,Irange,fbase)
% Fast image plotting, with scaling and colormap
% Adding image writing
I = double(I);
I = I - mean(I(:));
I = I / sqrt(mean(I(:).^2));
figure(1)
clf
imagesc(I)
axis equal off
colormap(violetFire(256))
set(gca,'position',[0 0 1 1])
if nargin == 1
    caxis([-2 2])  % default +/- two SDs of color?
elseif isempty(Irange)
    Irange = [min(I(:)) max(I(:))];
else
    caxis(Irange) 
end
if nargin == 3 || nargout > 0
    I = (I - Irange(1))/(Irange(2)-Irange(1));
    I(I<0) = 0;
    I(I>1) = 1;
end
if nargin == 3
    fname = ['gray_' fbase '_' num2str(Irange(1)) '_to_' ...
        num2str(Irange(2)) '.png'];
    imwrite(round(255*I)+1,gray(256),fname,'png');
    fname = [fbase '_' num2str(Irange(1)) '_to_' ...
        num2str(Irange(2)) '_color.png'];
    imwrite(round(255*I)+1,violetFire(256),fname,'png');
end

end