function res = ifftc_yj(x,dim)
%res = ifftc(x,dim)
res = fftshift(ifft(ifftshift(x,dim),[],dim),dim);
