%% initialization


foldername=[];filename=[];%category name and image name
Ensemble=[];%aggragate trainvector
Imagesize=[320 240]*2;%downsize image if larger than this size

%Imagesize=[640 480];%downsize image if larger than this size
%Graylevel=16;%gray level for cocurrence matrix
%Grid=[30 40];%patch size
offset=[0];%index offset of each category 
%texturedim=12;%dimension of texture feature in each patch
feature_dim=128;%dimension of feature in each patch
Ensemblesize=8e+5;
codewordnum=400;%size of codebook
% Num_patch=[10 10;
%            8 8;
%            15 15];%number of patches in an image
subpatchsize=48;%9,36,48;
MaxImagePatch=3;
%quantizelevel=2.^[1:8];
%quantizelevel=[32,128];
quantizelevel=[128];
stdskip=5;
filetype='.sift';%klentropy  sift  sifte siftentropy
% sift sifte sifta
mytree=[];