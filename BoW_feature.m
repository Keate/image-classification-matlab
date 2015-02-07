function BoW_feature(root,dataset)
initialization
Pathname=strcat(root,dataset,'\');%'D:\databases\101_ObjectCategories\'
Category=dir(strcat(Pathname,'*.*'));%read the database
CategoryNum=size(Category,1);%number of categories


%% feature extraction
for Ii=1:CategoryNum,

    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')), 

        foldername=Category(Ii).name;
        Image=dir(strcat(Pathname,foldername,'\*.jpg'));

        imgnum=1:length(Image);

        for k=imgnum,

    

               %entropyfile=strcat(Pathname,foldername,'\',Image(k).name,filetype);
                entropyfile=strcat(Pathname,foldername,'\',Image(k).name);dotpos=find(entropyfile=='.');
               entropyfile=strcat(entropyfile(1:dotpos(end)-1),'_sift.mat');%feature file is *.mat des, dsift, daisy
               %fr = fopen(entropyfile, 'r');
               %fr = -1;
               %if fr == -1,
                    I = imread(strcat(Pathname,foldername,'\',Image(k).name));Image(Ii).size=size(I);I=double(I);
                    fw = fopen(entropyfile, 'w');

                    if strcmp(filetype,'.sifte'),
                        descriptors=GenerateFastSiftDescriptors( I, 8, 16, 1 );
                    end
                    if strcmp(filetype,'.sift'),
                        descriptors=GenerateSiftDescriptors( I, 8, 16 ,1);descriptors = sp_normalize_sift(double(descriptors))*340;
                    end

                    %fwrite(fw, descriptors, 'double'); 
                    descriptors=single(descriptors);
                    save(entropyfile, 'descriptors', '-v7'); %v7 to compress the file 
                    fclose(fw);
               %else
               %     fclose(fr);
               %end

        end

        display(strcat(foldername,' feature extracted'))

    end


end


