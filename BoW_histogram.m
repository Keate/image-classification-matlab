function BoW_histogram(root,dataset)
initialization
%% obtain histograms for each image
%C=load('codebook_Caltech101.txt');
C=load('codebook_single.txt');
Pathname=strcat(root,dataset,'\');%'D:\databases\101_ObjectCategories\'
Category=dir(strcat(Pathname,'*.*'));%read the database
CategoryNum=size(Category,1);
%Itest=imread('D:\databases\NTU scene\Tiny NTU\flower\204902928_af39abcad7.jpg');
load EnsembleMeanVar 
Hist=[];
offset=[0];
for Ii=1:CategoryNum,
    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')), 
        foldername=Category(Ii).name;
        Image=dir(strcat(Pathname,foldername,'\*.jpg'));
        imgnum=0;
        for k=1:length(Image),
            enquiry=strcat(Pathname,foldername,'\',Image(k).name);
            Itest = imread(enquiry);Image(Ii).size=size(Itest);I=double(Itest);
            %Num_patch(1,1)=floor(size(I,1)/Grid(1));Num_patch(1,2)=floor(size(I,2)/Grid(2));Num_patch(2,:)=ceil(Num_patch(1,:)/2);Num_patch(3,:)=floor(Num_patch(1,:)*2);
            f=[];
            for kk=1:1, %multiresolution
                f=[f;ImagePatchFeatureHist(Itest,feature_dim,Num_patch(kk,:),Imagesize)];
                %f=[f;ImagePatchFeature(Itest,feature_dim,texturedim,Graylevel,Num_patch(kk,:),Imagesize)];
            end
            
            %f=(f-repmat(EnsembleMean,size(f,1),1))./repmat(EnsembleVar+eps,size(f,1),1);
            %fMean=mean(f);fVar=std(f); f=(f-repmat(fMean,size(f,1),1))./repmat(fVar+eps,size(f,1),1);
            
            %f=f.*repmat([1,ones(1,feature_dim-1)],size(f,1),1);
            index=zeros(1,size(f,1));
            for i=1:size(f,1),
                dis=sum((repmat(f(i,:),codewordnum,1)-C).^2,2);
                pos=find(dis==min(dis));
                index(i)=pos(1);
            end
            hist_test=hist(index,[1:codewordnum]);hist_test=hist_test/sum(hist_test);
            %bar(hist_test)
            imgnum=imgnum+1;
            Hist(offset(end)+imgnum,:)=hist_test;%histogram for every image
        end
        offset=[offset,offset(end)+imgnum];
    end
    display(strcat(foldername,' histogram calculated'))
end

totalimgnum=offset(end);
if strfind(lower(dataset),'test'),
    save Hist_test Hist CategoryNum totalimgnum
end
if strfind(lower(dataset),'train'),
    save Hist_train Hist CategoryNum totalimgnum
end