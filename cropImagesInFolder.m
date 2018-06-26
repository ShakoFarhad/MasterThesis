if(exist('imageFolder', 'var') && ischar(imageFolder) ~= 0)
    imageFolder = uigetdir(imageFolder, 'Select Image Folder');
else
    imageFolder = uigetdir('C:\', 'Select Image Folder');
end
addpath(imageFolder);
imageDir = dir(imageFolder);
fieldNames = fieldnames(imageDir);
images = rmfield(imageDir, fieldNames(2:end));
images = struct2cell(images);
images(:,1:2) = [];
images(end) = [];

figure
img = imread(images{end});
imagesc(img)
warning('off', 'Images:initSize:adjustingMag');
uiwait(helpdlg('Click, hold and drag to draw out a rectangle over the image.'));
rect = getrect;
x1 = round(rect(1));
x2 = round(x1 + rect(3));
y1 = round(rect(2));
y2 = round(y1 + rect(4));
if(x1 < 1)
    x1 = 1;
elseif(x1 > size(img,2))
    x1 = size(img,2);
end
if(y1 < 1)
    y1 = 1;
elseif(y1 > size(img,1))
    y1 = size(img, 1);
end
if(x2 < 1)
    x2 = 1;
elseif(x2 > size(img,2))
    x2 = size(img,2);
end
if(y2 < 1)
    y2 = 1;
elseif(y2 > size(img,1))
    y2 = size(img, 1);
end
close all
mkdir(imageFolder, 'croppedImages')

for imageNumber=1:size(images,2)
    img = imread(images{imageNumber});
    cropImg = img(y1:y2, x1:x2, :);
    %baseFileName = sprintf(images{imageNumber}, imageNumber);
    fullFileName = fullfile([imageFolder '\croppedImages'], images{imageNumber});
    imwrite(cropImg, fullFileName);
end