function [ImageName,ImageNum] = GetImageName(root, dataset)


Pathname=strcat(root,dataset,'\');
Category=dir(strcat(Pathname));
CategoryNum=size(Category,1);
ImageNum=0;
for Ii=1:CategoryNum,
    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')),
        foldername=Category(Ii).name;
        Image1=dir(strcat(Pathname,foldername,'\*.jpg'));
        Image2=dir(strcat(Pathname,foldername,'\*.png'));
        Image3=dir(strcat(Pathname,foldername,'\*.bmp'));
        Image=[Image1;Image2;Image3];
        ImageNum=ImageNum+length(Image);
    end
end


ImageName=cell(ImageNum,1);
imagenum=0;
for Ii=1:CategoryNum,
    if (Category(Ii).isdir==1 && ~strcmp(Category(Ii).name,'.') && ~strcmp(Category(Ii).name,'..')),
        foldername=Category(Ii).name;
        Image1=dir(strcat(Pathname,foldername,'\*.jpg'));
        Image2=dir(strcat(Pathname,foldername,'\*.png'));
        Image3=dir(strcat(Pathname,foldername,'\*.bmp'));
        Image=[Image1;Image2;Image3];
        for k=1:length(Image),
            imagenum=imagenum+1;
            ImageName{imagenum}=strcat(Pathname,foldername,'\',Image(k).name);
        end
    end
end