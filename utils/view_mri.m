function status = view_mri(varargin)
%%%% Yang 6/25/2021: show 3d matlab variables 

a=varargin{1};
%nii.img=rot90(fliplr(a),3);
nii.img=rot90(a,3);
nii.hdr=struct();
sz=size(nii.img);
if(length(sz)<3)
    sz(3)=1;   
end
nii.hdr.dime.dim=[3 sz 1 1 1];
nii.hdr.dime.datatype=0;
nii.hdr.dime.pixdim=[1 1 1 1 1 1];
nii.hdr.hist.originator=[sz/2 1];
view_nii(nii);
 
% nii=load_nii(a);
% status=view_nii(nii);