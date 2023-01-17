clc;clear;close all;
addpath /home/fs0/yj/wave_TOF/;
addpath /home/fs0/yj/wave_TOF/utils/;
addpath /home/fs0/yj/wave_TOF/utils/NIfTI_20140122/;
set(0,'DefaultFigurewindowStyle','docked')

%%
data_path='/vols/Scratch/hvl536/Data_TOF_Wave/MID00803';
MID='MID00803';
load(fullfile(data_path,'raw_data','CS_R4_mask.mat'));


fprintf('job start:CS_R4 \n')
tic

for ii=1:4
    
    if false%isfile(fullfile(data_path,[MID '_Result_CS_R8_Slab' num2str(ii) '.mat']))
        
        load (fullfile(data_path,'results',[MID '_Result_CS_R4_Slab' num2str(ii) '.mat']));
    else
        load(fullfile(data_path,'raw_data',[MID '_Raw_Slab' num2str(ii) '.mat']))
        [sx, sy, sz, sc] = size(sens_svd);
        num_chan=sc;
        psf_len=sx*6;
        m3d = repmat(permute(m2d_vd, [3,1,2]), [psf_len,1,1]);
        M3d = repmat(m3d, [1,1,1,num_chan]);
        
        % wave_encoded k-space and corresponding psf
        kdata_svd=mrir_zeropad(kdata_svd,[sx*6-size(kdata_svd,1) 0 0 0],'pre');
        psf_use  =  padarray(psf_use,[sx*6-size(psf_use,1) 0 0],1,'pre');
        data_Wave=M3d.*kdata_svd;%
        PSF_shift_Wave = repmat(psf_use, [1,1,1,num_chan]);
        
        % wave_decoded k-space and corresponding psf
        Img_KxYZ = ifftc(ifftc(kdata_svd,2),3);
        Img_KxYZ = Img_KxYZ .* repmat(conj(psf_use), [1,1,1,size(Img_KxYZ,4)]);
        kdata_dec= fftc(fftc(Img_KxYZ,2),3);
        psf_dec=ones(size(psf_use));
        data_noWave=M3d.*kdata_dec;%
        PSF_shift_noWave = repmat(psf_dec, [1,1,1,num_chan]);
        
        
        % check if wavelet toolbox exists
        v = ver;
        has_wavelet = any(strcmp(cellstr(char(v.Name)), 'Wavelet Toolbox'));
        
        
        %--------------------------------------------------------------------------
        % fnlCG recon for CS-Wave
        %--------------------------------------------------------------------------
        Nwav = [sx,sy,72];        % zero padded mtx size for wavelet transform
        
        param = init;
        
        % 3d-wavelet parameters:
        param.wav_scale = 5;
        param.wav_type = 'db1';
        param.wav_mode = 'per';
        
        
        param.TV = TV3D;
        param.M3d = M3d;
        
        param.TVWeight = 3e-4;                      % 3D TV penalty
        param.WavWeight = has_wavelet * 3e-4;       % 3D Wavelet penalty: set to 0 if wavelet toolbox not present
        
        param.num_chan = num_chan;
        param.crop_size = sx;
        param.pad_size = (psf_len - sx) / 2;
        param.wav_pad = (Nwav - [sx sy sz]) / 2;
        
        param.Receive = sens_svd;
        param.Ct = conj(param.Receive);
        
        param.fft_scale = sqrt(psf_len * sy * sz);
        
        param.data = data_Wave*param.fft_scale;
        param.PSF_shift = PSF_shift_Wave;
        param.PSFt_shift = conj(PSF_shift_Wave);
        
        param.Itnlim = 150;
        param.tol = 1e-3;
        param.pNorm = 1;
        
        
        % run cs_wave
        img_WaveCS = zeros(sx,sy,sz);
        tic
        img_WaveCS = fnlCg_ics_wave_yj(img_WaveCS, param);
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
        
        
        % run cs_STD
        param.data = data_noWave*param.fft_scale;
        param.PSF_shift = PSF_shift_noWave;
        param.PSFt_shift = conj(PSF_shift_noWave);
        img_CS = zeros(sx,sy,sz);
        tic
        img_CS = fnlCg_ics_wave_yj(img_CS, param);
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))

    end
    

    
end

toc
fprintf('job end \n')

