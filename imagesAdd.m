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
addpath(imageFolder);
image = dir(imageFolder);
fieldNames = fieldnames(image);
images = rmfield(image, fieldNames(2:end));
images = struct2cell(images);
images(:,1:2) = [];
images(end) = [];
amountOfImagesBeforeRemoval = size(images,2); %size before removal.
prompt = {'Skip every X image:'};
dlg_title = 'Creating Masks';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if(isempty(answer))
    skipImage = 0;
else
    skipImage = str2double(cell2mat(answer(1)));
end

if(skipImage ~= 0 && skipImage > 1)
    images(:,1:skipImage:end) = [];
end
amountOfImagesAfterRemoval = size(images,2); %Size after removal
if(amountOfImagesBeforeRemoval ~= amountOfImagesAfterRemoval)
    fprintf('Amount of images before and after; %.f, %.f. \n', amountOfImagesBeforeRemoval, amountOfImagesAfterRemoval);
    fprintf('%.f images will not be included in the calculations.\n', amountOfImagesBeforeRemoval - amountOfImagesAfterRemoval);
end

%Initialization of variables
tempImage1 = imread(images{1});
tempImage2 = imread(images{2});
if(size(tempImage1,3) == 3)
    bilde1 = im2double(rgb2gray(tempImage1));
else
    bilde1 = im2double(tempImage1);
end
if(size(tempImage1,3) == 3)
    bilde2 = im2double(rgb2gray(tempImage2));
else
    bilde2 = im2double(tempImage2);
end
imageAdd = imadd(bilde1,bilde2);

for imageNumber=2:size(images,2)-1
    tempImage1 = imread(images{imageNumber});
    if(size(tempImage1,3) == 3)
        bilde = im2double(rgb2gray(tempImage1));
    else
        bilde = im2double(tempImage1);
    end
    imageAdd = imadd(bilde, imageAdd);
end

retry = 1;
fontSize = 30;
startValues = 1;
while(retry == 1)
    if(startValues == 1)
        sigma = 0.5;
        threshold = 1;
        alphaUpper = 1.3;
        alphaLower = 1.3;
        startValues = 0;
    else
        prompt = {'Blur Strength:','Black and White Threshold:', 'Percentage upper mask cutoff', 'Percentage lower mask cutoff'};
        dlg_title = 'Create Mask From Images';
        num_lines = 1;
        defaultans = {num2str(sigma),num2str(threshold), num2str(alphaUpper),num2str(alphaLower)};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if(~isempty(answer))
            %This is for the gauss filter.
            sigma = str2double(cell2mat(answer(1)));
            %Threshold for making the image binary
            threshold = str2double(cell2mat(answer(2)));
            %Alpha is the percentage mask cutoff. It determines how much of the mask
            %shall be discarded.
            alphaUpper = str2double(cell2mat(answer(3)));
            alphaLower = str2double(cell2mat(answer(4)));
        end
    end

    %Making the image black and white, and inversing the colors
    imageAddPostProssessed = min(imageAdd,mean(mean(imageAdd))*1.5);
    imageAddPostProssessed = imageAddPostProssessed/mean(mean(imageAddPostProssessed));
    if(sigma ~= 0)
        imageAddPostProssessed = imgaussfilt(imageAddPostProssessed,sigma);
    end
    
    normalized = imageAddPostProssessed;
    normalized(imageAddPostProssessed < threshold) = 0;
    normalized(imageAddPostProssessed >= threshold) = 1;

    maskDistances = bwdist(imcomplement(normalized));
    [maxDistances, y] = max(maskDistances);
    x = 1:size(y,2);
    for i=x
        normalized(y(i)+round(maxDistances(i)*alphaLower):size(normalized,1),i) = 0;
        normalized(1:y(i)-round(maxDistances(i)*alphaUpper) ,i) = 0;
    end
    complementNormalized = imcomplement(normalized);
    
    f = figure();
    subplot(1,2,1);
    imagesc(imageAddPostProssessed)
    title('Added Image Intensities','FontSize',fontSize)
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label
    set(gca, 'FontSize', fontSize);
    subplot(1,2,2);
    imagesc(complementNormalized)
    colormap('gray')
    title('Added Image Intensities Vein Mask','FontSize',fontSize)
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label
    set(gca, 'FontSize', fontSize);
    aspectRatio = min(size(normalized,2))/min(size(normalized,1));
    set(f, 'Position', [-100, 150, max(min(size(normalized,2)*3.5,700*2.5*aspectRatio),600), max(min(size(normalized,1),700),350)])
    
    maxValue = max(max(imageAddPostProssessed));
    minValue = min(min(imageAddPostProssessed));
    answer = questdlg(sprintf('Choose Threshold Between: [%.2f, %.2f]', minValue, maxValue), ...
        'Is the mask good?');
    % Handle response
    switch answer
        case 'Yes'
            retry = 0;
        case 'No'
            retry = 1;
            close all;
        case 'Cancel'
            retry = 2;
            close all
    end
end
if(retry == 0)
    if(exist('maskFolder', 'var') && ischar(maskFolder) ~= 0)
        maskFolder = uigetdir(maskFolder, 'Select Mask Folder');
    else
        maskFolder = uigetdir(imageFolder, 'Select Mask Folder');
    end
    
    imageName = [maskFolder, '\MaskWithAxis'];
    print(f,imageName,'-dpng');

    imageName = [maskFolder, '\Mask.png'];
    imwrite(complementNormalized,imageName);
    
    imageName = [maskFolder, '\InvertedMask.png'];
    imwrite(normalized,imageName);
    close all
end