function TrainTestPartition(Pathname)

%Pathname='D:\databases\Caltech 101\';
Category=dir(strcat(Pathname,'*.*'));%read the database
CategoryNum=size(Category,1);%number of categories
%% feature extraction
for Ii=1:CategoryNum,
    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')), 
        foldername=Category(Ii).name;
        Image=dir(strcat(Pathname,foldername,'\*.jpg'));
        index=randperm(length(Image));
        mkdir(strcat(Pathname,'train\',foldername,'\'));
        for k=1:floor(length(Image)/2),
            movefile(strcat(Pathname,foldername,'\',Image(index(k)).name),strcat(Pathname,'train\',foldername,'\',Image(index(k)).name),'f');
        end
        mkdir(strcat(Pathname,'test\',foldername,'\'));
        for k=floor(length(Image)/2)+1:length(Image),
            movefile(strcat(Pathname,foldername,'\',Image(index(k)).name),strcat(Pathname,'test\',foldername,'\',Image(index(k)).name),'f');
        end
        rmdir(strcat(Pathname,foldername,'\'),'s')
        display(strcat(foldername,' partitioned'))
    end
end