clc;clear;close all;
addpath /home/fs0/yj/wave_TOF/utils/;
addpath /home/fs0/yj/wave_TOF/utils/NIfTI_20140122/;
set(0,'DefaultFigurewindowStyle','docked')

%%
data_path='/vols/Scratch/hvl536/Data_TOF_Wave/MID00803';
MID='MID00803';

fprintf('job start:CAIPI_R4\n')
tic
for ii=1:4
    
    if false%isfile(fullfile(data_path,[MID '_Result_CS_R4_Slab' num2str(ii) '.mat']))
        
        load(fullfile(data_path,'results',[MID '_Result_Wave_R4_Slab' num2str(ii) '.mat']));
    else
        load(fullfile(data_path,'raw_data',[MID '_Raw_Slab' num2str(ii) '.mat']))
        [sx, sy, sz, sc] = size(sens_svd);
        num_chan=sc;
        psf_len=sx*6;
        
        
        %mask
        Ry=4;Rz=1;
        mask_yz_tmp = zeros(sy, sz); % sampling mask
        mask_yz_tmp(1:Ry:end, 1:Rz:end) = 1;
        for i=1:sz
            shift_length=2*(mode(i,4)-1);
            mask_yz_tmp(:,i)=circshift(mask_yz_tmp(:,i),shift_length);
        end   
        mask_yz=mask_yz_tmp(1:256,:);
%         figure, imshow(mask_yz'), title('sampling mask'); % 2D CAIPI pattern
        mask = repmat(reshape(mask_yz, [1, sy, sz]), [psf_len, 1, 1,num_chan]);
        
        
        
        % wave_encoded k-space and corresponding psf
        kdata_svd=mrir_zeropad(kdata_svd,[sx*6-size(kdata_svd,1) 0 0 0],'pre');
        psf_use  =  padarray(psf_use,[sx*6-size(psf_use,1) 0 0],1,'pre');
        
        
        % wave_decoded k-space and corresponding psf
        Img_KxYZ = ifftc(ifftc(kdata_svd,2),3);
        Img_KxYZ = Img_KxYZ .* repmat(conj(psf_use), [1,1,1,size(Img_KxYZ,4)]);
        kdata_dec= fftc(fftc(Img_KxYZ,2),3);
        psf_dec=ones(size(psf_use));
        
        
        % parameter
        param = [];
        param.psf_len = psf_len;                             % psf length
        param.img_len = sx;                                    % image length (without oversampling)
        param.pad_size = (param.psf_len - param.img_len)/2;             % (psf length - image length) / 2
        param.num_chan = num_chan;
        param.Nx=sx;
        param.Ny=sy;
        param.Nz=sz;
        param.psfs = repmat(psf_use, [1,1,1,param.num_chan]);
        param.sens = sens_svd;
        param.m3d = mask;


        lsqr_iter = 200;
        lsqr_tol = 1e-3;
        
        
        % run wave_caipi
        kdata_svd = kdata_svd.*mask;
        tic;
        img=lsqr(@apply_wave_yj, kdata_svd(:), lsqr_tol, lsqr_iter, [], [], [], param);
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
        img_wave=reshape(img,[sx sy sz]);
        
        % run std_caipi
        param.psfs = repmat(psf_dec, [1,1,1,param.num_chan]);
        kdata_dec = kdata_dec.*mask;
        tic;
        img=lsqr(@apply_wave_yj, kdata_dec(:), lsqr_tol, lsqr_iter, [], [], [], param);
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
        img_std=reshape(img,[sx sy sz]);

    end
    
    
end
toc
fprintf('job end \n')