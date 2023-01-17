function [ scale_axis ] = genCaxis( img, std_factor )
%GENCAXIS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    std_factor = 4;
end

if norm(imag(img(:))) > 0
    % complex => convert to magnitude image
    img_avg = mean( abs(img(img~=0)) );
    img_std = std( abs(img(img~=0)) );
    
    scale_axis = [0, img_avg + std_factor*img_std];
else
    % real
    if norm( img(img<0) ) > 0
        % has negative values
        img_avg = mean( img(img~=0) );
        img_posStd = std( img(img>0) );
        img_negStd = std( img(img<0) );

        scale_axis = [img_avg + std_factor*img_negStd, img_avg + std_factor*img_posStd];
    else
        % magnitude image
        img_avg = mean( img(img~=0) );
        img_std = std( img(img~=0) );

        scale_axis = [0, img_avg + std_factor*img_std];
    end
end

end

