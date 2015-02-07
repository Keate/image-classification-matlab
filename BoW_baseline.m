clear,clc,close all,fclose all


%cd '.\gbvs'
%gbvs_install;
%cd '..'



cd '.\toolbox'
vl_setup;
cd '..'
initialization

%% choose the database path
%root ='\\eee-mtl-dellykh\YKH_PhDStudents\Zhang Dajiang\SingaporeLandmark';
%root ='D:\databases\scene_categories';
%root ='D:\tagSuggestionNtuWideS\NTU scene';
%root ='D:\databases\Caltech 8';
%root = 'D:\leaf';
%root = 'D:\oxbuild_images';
%root = 'D:\PFIs';
%root = 'I:\study\mypapers\under work\oxbuild_images';
%root = 'D:\under work\object localization\leaf dataset\leaf3\';
root = 'D:\under work\object localization\leaf dataset\leaf5\';
%root = '/media/My Passport/study/mypapers/under\ work/BoWBaseline_saliency';
%% compute the feature vectors for codebook generation 

%tic,dataset='\Train';BoW_feature(root,dataset),toc
%tic,dataset='\Test';BoW_feattture(root,dataset),toc

tic,dataset='\Train';BoW_feature2(root,dataset),toc
tic,dataset='\Test';BoW_feature2(root,dataset),toc
% do clustering to generate codebook
tic,
MyKmeans(root);
toc

%% compute the histograms of training images
dataset='\Train';tic,BoW_GlobalSIFT2(root,dataset),toc

%% compute the histograms of test images
dataset='\Test';tic,BoW_GlobalSIFT2(root,dataset),toc

%% train the classifiers

%dataset='\Train';tic,fasttrain(root,dataset),toc


%% classify test images

%dataset='\Test';tic,fastclassify(root,dataset),classificationtime=toc

% dataset='\Test';
% %tic, Geometric_verify(root,dataset,1,0,1), toc
tic, Geometric_verify0(root,dataset,0,0,0,1), toc  %root,dataset,b_siftmatch, b_ransac, b_selkdtree, topN 
%  tic, oxford_query(root,dataset,0,0,0,1), toc
% %tic, Geometric_verify0(root,dataset,0,0,0,1), toc %without GV

