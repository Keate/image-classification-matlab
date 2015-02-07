function BoW_feature2(root,dataset)
initialization
Pathname=strcat(root,dataset,'\');%'D:\databases\101_ObjectCategories\'
Category=dir(strcat(Pathname));%read the database
CategoryNum=size(Category,1);%number of categories

RandCategoryNum=randperm(CategoryNum);
%% feature extraction
for Ii=RandCategoryNum %1:CategoryNum,

    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')), 

        foldername=Category(Ii).name;
        Image1=dir(strcat(Pathname,foldername,'\*.jpg'));
        Image2=dir(strcat(Pathname,foldername,'\*.png'));
        Image3=dir(strcat(Pathname,foldername,'\*.bmp'));
        Image=[Image1;Image2;Image3];

        %imgnum=1:length(Image);
        imgnum=randperm(length(Image));

        for k=imgnum,
               
               entropyfile=strcat(Pathname,foldername,'\',Image(k).name);dotpos=find(entropyfile=='.');
               entropyfile=strcat(entropyfile(1:dotpos(end)-1),'_sift.mat');%feature file is *.mat des, dsift, daisy
               if ~exist(entropyfile),
                    I = imread(strcat(Pathname,foldername,'\',Image(k).name));Image(Ii).size=size(I);I=double(I);
                    if strcmp(filetype,'.sift'),
                        %[descriptors,des_weight]=GenerateSiftDescriptors( I, 8, 16, 1);
%                         descriptors=bsxfun(@rdivide,descriptors,sum(descriptors.^2,2).^(1/2));
%                         descriptors=bsxfun(@rdivide,descriptors,max(descriptors')')*255;
                        [descriptors]=GenerateFastSiftDescriptors2( I, 8, 16, 1 );
                        %[descriptors, des_weight]=GenerateDaisyDescriptors( I, 8, 16, 1 , 1);
                        descriptors=single(descriptors);
                        %save(entropyfile, 'descriptors', 'des_weight', '-v7'); %v7 to compress the file      
                        save(entropyfile, 'descriptors', '-v7'); %v7 to compress the file 
                    end
               end

        end

        display(strcat(foldername,' feature extracted'))

    end


end


