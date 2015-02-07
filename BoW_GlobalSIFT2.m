function BoW_GlobalSIFT2(root,dataset)
initialization
%% obtain histograms for each image

%C=load('codebook_GlobalSIFT.txt');
load mytree
%load mytree1
Pathname=strcat(root,dataset,'\');
Category=dir(strcat(Pathname));
CategoryNum=size(Category,1);

%load EnsembleMeanVar

if ~isempty(strfind(lower(dataset),'test')),
    %Hist=zeros(200,10^4,'single');
    Hist=zeros(200,10101,'single');
end
if ~isempty(strfind(lower(dataset),'train')),
    %Hist=zeros(200,10^4,'single');
    Hist=zeros(200,10101,'single');
end

%Hist=[];
offset=[0];
ImageName={};
for Ii=1:CategoryNum,
    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')),
        foldername=Category(Ii).name;
        Image1=dir(strcat(Pathname,foldername,'\*.jpg'));
        Image2=dir(strcat(Pathname,foldername,'\*.png'));
        Image3=dir(strcat(Pathname,foldername,'\*.bmp'));
        Image=[Image1;Image2;Image3];
        imgnum=0;
        for k=1:length(Image),
            entropyfile=strcat(Pathname,foldername,'\',Image(k).name);dotpos=find(entropyfile=='.');
            entropyfile=strcat(entropyfile(1:dotpos(end)-1),'_sift.mat');%feature file is *.mat
            load(entropyfile);%descriptors=descriptors*310;



            descriptorspath= vl_hikmeanspush(mytree,uint8(descriptors'));
            %descriptorspath= vl_ikmeanspush(uint8(descriptors'),mytree);

            hist_test = vl_hikmeanshist(mytree,descriptorspath)';
            %hist_test=hist(descriptorspath,[1:size(mytree,2)]);
            %hist_test = single(vl_hikmeanshist2(mytree,descriptorspath,des_weight));


            imgnum=imgnum+1;
            if size(Hist,1)<offset(end)+imgnum,
                Hist=[Hist;hist_test];
            else
                Hist(offset(end)+imgnum,:)=hist_test;%histogram for every image
            end
        end
        offset=[offset,offset(end)+imgnum];
    end
    display(strcat(foldername,' histogram calculated'))
end
yapp = zeros(offset(end),1);
for i=1:CategoryNum-2,
    pos=offset(i+1);
    yapp(offset(i)+1:pos) = i;
end
categoryindex = yapp;

totalimgnum=offset(end);
if ~isempty(strfind(lower(dataset),'test')),
    testoffset=offset;
    save Hist_test Hist CategoryNum totalimgnum testoffset categoryindex ImageName
end
if ~isempty(strfind(lower(dataset),'train')),
    trainoffset=offset;
    save Hist_train Hist CategoryNum totalimgnum trainoffset categoryindex ImageName
end
