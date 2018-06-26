%This needs to be added to run the Hydrolav PIV code.
p = genpath('HydrolabPIV\src');
addpath(p);
javaaddpath(['HydrolabPIV\src' filesep 'measures']);
javaaddpath(['HydrolabPIV\src' filesep 'interp']);

if(exist('imageFolder', 'var') && ischar(imageFolder) ~= 0)
    try
        imageFolder = uigetdir(imageFolder, 'Select Image Folder');
    catch
        disp('No image folder selected')
        return
    end
else
    try
        imageFolder = uigetdir(pwd, 'Select Image Folder');
    catch
        disp('No image folder selected')
        return
    end
end
if(ischar(imageFolder) == 0)
    disp('No image folder selected')
    return
end
addpath(imageFolder);
imageDir = dir(imageFolder);
fieldNames = fieldnames(imageDir);
images = rmfield(imageDir, fieldNames(2:end));
images = struct2cell(images);
images(:,1:2) = [];
images(end) = [];
amountOfImagesBeforeRemoval = size(images,2); %size before removal.

prompt = {'Skip every X image:', 'Use subset of images: Startpoint', 'Use subset of images: Endpoint'};
dlg_title = 'Ensemble PIV with Distorted Pass';
num_lines = 1;
defaultans = {'0', '1', num2str(amountOfImagesBeforeRemoval)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

skipImage = str2double(cell2mat(answer(1)));
startPoint = str2double(cell2mat(answer(2)));
endPoint = str2double(cell2mat(answer(3)));

if(skipImage ~= 0 && skipImage > 1)
    images(:,1:skipImage:end) = [];
end
amountOfImagesAfterRemoval = size(images,2); %Size after removal
if(startPoint >= 1 && endPoint <= amountOfImagesBeforeRemoval && endPoint <= amountOfImagesAfterRemoval)
    images = images(:,startPoint:endPoint);
end
amountOfImagesAfterRemoval = size(images,2); %Size after removal
if(amountOfImagesBeforeRemoval ~= amountOfImagesAfterRemoval)
    fprintf('Amount of images before and after; %.f, %.f. \n', amountOfImagesBeforeRemoval, amountOfImagesAfterRemoval);
    fprintf('%.f images will not be included in the calculations.\n', amountOfImagesBeforeRemoval - amountOfImagesAfterRemoval);
end

%Setting up the loading dialog window
waitBar = waitbar(0,'Please wait... PIV 0.0% finished.','Name','Ensemble PIV with Distorted Pass');

overlap = 0.75;
x = 240;
y = 60;
x1 = round(x/3);
y1 = round(y/3);

opt1 = setpivopt('range',[-x1 x1 -y1 y1],'subwindow',x,y,overlap,'savepeaks',true);
imagesUsedInEnsemblePiv = 0;
for imageNumber=1:size(images,2)-1
    try
        loadingPercentText = ['Please wait... Ensemble PIV ', sprintf('%.1f',imageNumber*100.0 / (size(images,2)*2)), '% finished.'];
        waitbar(imageNumber / (size(images,2)*2), waitBar, loadingPercentText)
    catch
        disp("Loading bar cancelled. Exiting PIV.");
        imagesUsedInEnsemblePiv = imageNumber;
        break
    end

    tempImage1 = imread(images{imageNumber});
    tempImage2 = imread(images{imageNumber + 1});
    if(size(tempImage1,3) == 3)
        im1 = im2double(rgb2gray(tempImage1));
    else
        im1 = im2double(tempImage1(:,:,1));
    end
    if(size(tempImage1,3) == 3)
        im2 = im2double(rgb2gray(tempImage2));
    else
        im2 = im2double(tempImage2(:,:,1));
    end
    if(imageNumber==1)
        piv1 = normalpass([],im1,[],im2,[],opt1);
    else
        piv1 = normalpass([],im1,[],im2,[],piv1);
    end
end


if(imagesUsedInEnsemblePiv < size(images,2)-1)
    %Setting up the loading dialog window
    waitBar = waitbar(0,['Please wait... PIV', sprintf('%.1f',imageNumber*100.0 / (size(images,2)*2)), ' finished.'],'Name','Ensemble PIV with Distorted Pass');
end
opt2=setpivopt('range',[-x1 x1 -y1 y1],'subwindow',x,y,overlap,'savepeaks',true);
for imageNumber=1:size(images,2)-1
    try
        loadingPercentText = ['Please wait... Distorted Pass PIV ', sprintf('%.1f',(imagesUsedInEnsemblePiv + imageNumber)*100.0 / (imagesUsedInEnsemblePiv + size(images,2))), '% finished.'];
        waitbar((imagesUsedInEnsemblePiv + imageNumber) / (imagesUsedInEnsemblePiv + size(images,2)), waitBar, loadingPercentText)
    catch
        disp("Loading bar cancelled. Exiting PIV.");
        break
    end

    tempImage1 = imread(images{imageNumber});
    tempImage2 = imread(images{imageNumber + 1});
    if(size(tempImage1,3) == 3)
        im1 = im2double(rgb2gray(tempImage1));
    else
        im1 = im2double(tempImage1(:,:,1));
    end
    if(size(tempImage1,3) == 3)
        im2 = im2double(rgb2gray(tempImage2));
    else
        im2 = im2double(tempImage2(:,:,1));
    end

    if(imageNumber==1)    
        piv2 = distortedpass(piv1,im1,[],im2,[],opt2);
    else
        piv2 = distortedpass(piv1,im1,[],im2,[],piv2);
    end
end

try
    close(waitBar)
catch
    disp("No loading bar to close. Exiting PTV")
end

[U2,V2,x2,y2] = replaceoutliers(piv1);

fontSize = 30;
figure
load wind
cav = curl(x2,y2,U2,V2);%plotting vorticity
pcolor(x2,y2,cav); shading interp
hold on;
quiver(x2,y2,U2,V2,'y')
hold off
colormap winter
c = colorbar('southoutside','FontSize',fontSize);
c.Label.String = 'Velocity in Pixels per Frame';
set(gca, 'YDir', 'reverse');
set(gca, 'FontSize', fontSize);
titleText = ['PIV Ensemble Distorted Passes.'];
title(titleText,'FontSize',fontSize)
xlabel('Pixels','FontSize',fontSize) % x-axis label
ylabel('Pixels','FontSize',fontSize) % y-axis label