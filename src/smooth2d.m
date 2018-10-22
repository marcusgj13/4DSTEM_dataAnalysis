function [smoothImg cutoffRad]= smooth(img,resolutionCutoff)
    %2D version of smooth3D
    Rsize = size(img,1);
    Csize = size(img,2);
    Rcenter = round((Rsize+1)/2);
    Ccenter = round((Csize+1)/2);
    a=1:1:Rsize;
    b=1:1:Csize;
    [bb,aa]=meshgrid(b,a);
    sigma=(Rsize*resolutionCutoff)/(2*sqrt(2));
    kfilter=exp( -( ( ((sqrt((aa-Rcenter).^2+(bb-Ccenter).^2)).^2) ) ./ (2* sigma.^2) ));
    kfilter=kfilter/max(max(kfilter));
    kbinned = my_fft(img);  
   
    kbinned = kbinned.*kfilter;
    smoothImg = my_ifft(kbinned);
    
    [Y X] = ind2sub(size(img),find(kfilter<(exp(-1))));

    Y = Y-(size(img,2)/2);
    X = X-(size(img,2)/2);
    R = sqrt(Y.^2+X.^2);
    cutoffRad = ceil(min(abs(R)));
end