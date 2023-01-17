function [ res, tflag ] = apply_wave_yj( in, params, tflag )
% wave sense reconstruction using forward model
% Yang Ji, FMRIB, 2022
if strcmp(tflag,'transp')
    
    kspace = reshape(in, [params.psf_len,params.Ny,params.Nz,params.num_chan]);    
    kspace   = kspace .* params.m3d;% sampling mask;
    
    kspace    = ifftshift(ifft(fftshift(kspace,3), [], 3), 3).* sqrt(params.Nz);  %ifft along z dimension
    kspace     = ifftshift(ifft(fftshift(kspace ,2), [], 2), 2).* sqrt(params.Ny);  %ifft along y dimension
    img_coils    = ifftshift(ifft(fftshift(kspace .* conj(params.psfs),1),[],1),1).* sqrt(params.psf_len); %wave decoding
    img          = img_coils(end/2+1-params.img_len/2:end/2+params.img_len/2,:,:,:).* conj(params.sens);                                                                                               
    Res          = sum(img,4);
    res          = Res(:);
  
else
    
    % Forward wave SENSE operator
    img_coils = repmat(reshape(in, params.Nx,params.Ny,params.Nz), [1,1,1,params.num_chan]) .* params.sens;% 
    img_coils_oversamp = padarray(img_coils, [params.pad_size,0,0,0]);
    
    kspace      = fftshift(fft(ifftshift(img_coils_oversamp,1), [], 1), 1) .* params.psfs / sqrt(params.psf_len); %wave encoding
    kspace      = fftshift(fft(ifftshift(kspace ,2), [], 2), 2)/ sqrt(params.Ny);  % fft along y dimension
    kspace      = fftshift(fft(ifftshift(kspace,3), [], 3), 3)/ sqrt(params.Nz);  % fft along z dimension
    
    kspace    = kspace.*params.m3d; % sampling mask
    
    res           = kspace(:);
    
end

end
