function MIP = mip(array,dim,oversample)
%   mip(array, dim, oversample):
%   Calculates the maximum intensity projection (MIP) of array along the
%   dimension dim. In the case of complex data, it looks for the values
%   with the greatest magnitude. Oversample indicates the (factor of, in
%   case of <0) number of slices to remove, split over the outsides of
%   array in dimension dim (defaults to 0). Currently only works for a 
%   maximum of 6 dimensions.
%   matthijs.debuck@ndcn.ox.ac.uk
    if ~exist('oversample','var')
        oversample = 0;
    end
    
    if (oversample > 0) * (oversample < 1)
        oversample = int8(size(array,dim)*oversample);
    end
    
    oversample_top = oversample/2 - mod(oversample,2);
    oversample_bot = oversample/2;

    switch dim
        case 1
            MIP = squeeze(max(array(1+oversample_top:end-oversample_bot,:,:,:,:,:),[],dim));
        case 2
            MIP = squeeze(max(array(:,1+oversample_top:end-oversample_bot,:,:,:,:),[],dim));
        case 3
            MIP = squeeze(max(array(:,:,1+oversample_top:end-oversample_bot,:,:,:),[],dim));
        case 4
            MIP = squeeze(max(array(:,:,:,1+oversample_top:end-oversample_bot,:,:),[],dim));
        case 5
            MIP = squeeze(max(array(:,:,:,:,1+oversample_top:end-oversample_bot,:),[],dim));
        case 6
            MIP = squeeze(max(array(:,:,:,:,:,1+oversample_top:end-oversample_bot),[],dim));
    end
end

