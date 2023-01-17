function res = fftc_yj(x,dim)
% res = fftc(x,dim)
res = fftshift(fft(ifftshift(x,dim),[],dim),dim);
