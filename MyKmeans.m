function MyKmeans(root)

initialization

dataset='\Train';
Pathname=strcat(root,dataset,'\');%'D:\databases\101_ObjectCategories\'
Category=dir(strcat(Pathname));%read the database
CategoryNum=size(Category,1);%number of categories

Ensemble=zeros(Ensemblesize,feature_dim,'uint8');
%EnsembleWeight=zeros(Ensemblesize,1,'single');

C=[];
imagenum=0;
%% feature extraction
for Ii=1:CategoryNum,

    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')),

        foldername=Category(Ii).name;
        Image1=dir(strcat(Pathname,foldername,'\*.jpg'));
        Image2=dir(strcat(Pathname,foldername,'\*.png'));
        Image3=dir(strcat(Pathname,foldername,'\*.bmp'));
        Image=[Image1;Image2;Image3];

        for k=1:length(Image),
            entropyfile=strcat(Pathname,foldername,'\',Image(k).name);dotpos=find(entropyfile=='.');
            entropyfile=strcat(entropyfile(1:dotpos(end)-1),'_sift.mat');%feature file is *.mat
            load(entropyfile);%descriptors=descriptors*310;

               
            des_num_per_image=floor(size(Ensemble,1)/(CategoryNum-2)/length(Image));
            patchsize=min(des_num_per_image,size(descriptors,1));
            
            des_weight=ones(1,patchsize)/patchsize;
            [sorted_weight,patchindex]=sort(des_weight,'descend');


            patchindex=randperm(size(descriptors,1));
            descriptors=descriptors(patchindex(1:patchsize),:);
            
            Ensemble(imagenum+1:imagenum+patchsize,:)=uint8(descriptors);
            %EnsembleWeight(imagenum+1:imagenum+patchsize)=des_weight(patchindex(1:patchsize));
            


            imagenum=imagenum+patchsize;
        end
        display(strcat(foldername,' feature loaded'))
    end


end

%EnsembleWeight=EnsembleWeight(boolean(EnsembleWeight));
%Ensemble=Ensemble(1:length(EnsembleWeight),:);
Ensemble=Ensemble(boolean(sum(Ensemble,2)),:);
save Ensemble Ensemble
%save EnsembleWeight EnsembleWeight

load Ensemble
%load EnsembleWeight

% codevectorsize=size(Ensemble,1)
%
%
%     EnsembleMean=mean(Ensemble);
%     EnsembleMean=round(EnsembleMean);
%
% EnsembleStd=std(Ensemble);
% save EnsembleMeanVar EnsembleMean EnsembleStd
%  Ensemble=(Ensemble-repmat(EnsembleMean,size(Ensemble,1),1))./repmat(EnsembleStd,size(Ensemble,1),1);
%Ensemble=Ensemble-repmat(EnsembleMean,size(Ensemble,1),1)+max(EnsembleMean);
display('begin clustering')
tic

%[C0, A, ENERGY]=vl_kmeans(Ensemble', codewordnum, 'verbose', 'distance', 'l1', 'algorithm', 'elkan') ;C0=C0';
%C=double(C0);

% index=randperm(size(Ensemble,1));%random select centers
% index=index(1:11111);
% mytree=int32(Ensemble(index,:))';

[mytree,asgn] = vl_hikmeans(uint8(Ensemble'), 100, 100^2, 'method', 'elkan');
%[mytree,asgn] = vl_ikmeans(uint8(Ensemble'), 10^4, 'method', 'elkan');
%mytree = vl_myhikmeans(10, 4, EnsembleWeight);
save mytree mytree
%save mytree1 mytree1
display('end clustering, minutes:')
toc/60






% [C0, A, ENERGY]=vl_kmeans(Ensemble', codewordnum, 'verbose', 'distance', 'l1', 'algorithm', 'elkan') ;C0=C0';
% 
% C=double(C0);
% 
% display('end clustering, minutes:')
% 
% save codebook_GlobalSIFT.txt C -ascii