if(~exist('uniqueDataNameCounter', 'var'))
    %This needs to be added to run the Hydrolav PIV code at least once.
    p = genpath('HydrolabPIV\src');
    addpath(p);
    javaaddpath(['HydrolabPIV\src' filesep 'measures']);
    javaaddpath(['HydrolabPIV\src' filesep 'interp']);
    javaaddpath(['HydrolabPIV\src' filesep 'subpixel']);

    uniqueDataNameCounter = 1;
    settings = [];
else
    if(~exist(sprintf('settings_%d',uniqueDataNameCounter), 'var'))
        settings = [];
    else
        settings = eval(sprintf('settings_%d',uniqueDataNameCounter));
        particleLength = settings.particleLength;
        filter = settings.filter;
        noiseLength = settings.noiseLength;
    end
end

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

prompt = {'Skip every X image:', 'Use subset of images: Startpoint', 'Use subset of images: Endpoint','Particle Diameter in Pixels:', 'Image Filtering; Aggressive, Gaussian, Box, No:', 'Image Noise Length:'};
dlg_title = 'Particle Finder';
num_lines = 1;
if(isempty(settings))
    defaultans = {'0', '1', num2str(amountOfImagesBeforeRemoval), '4.0', 'Aggressive', '1.0'};
else
    defaultans = {'0', '1', num2str(amountOfImagesBeforeRemoval), num2str(particleLength), filter, num2str(noiseLength)};
end
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if(isempty(answer))
    skipImage = str2double(cell2mat(defaultans(1)));
    startPoint = str2double(cell2mat(defaultans(2)));
    endPoint = str2double(cell2mat(defaultans(3)));
    particleLength = str2double(cell2mat(defaultans(4))); 
    filter = defaultans(5); noiseLength = str2double(cell2mat(defaultans(6)));
else
    skipImage = str2double(cell2mat(answer(1)));
    startPoint = str2double(cell2mat(answer(2)));
    endPoint = str2double(cell2mat(answer(3)));
    particleLength = str2double(cell2mat(answer(4))); 
    filter = answer(5); noiseLength = str2double(cell2mat(answer(6)));
end

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

myGUI(images, imageFolder, particleLength, filter, noiseLength, uniqueDataNameCounter, settings)

function myGUI(images, imageFolder, particleLength, filter, noiseLength, uniqueDataNameCounter, settings)
    tempImage = imread(images{1});
    if(size(tempImage,3) == 3)
        firstImage = im2double(rgb2gray(tempImage));
    else
        firstImage = im2double(tempImage(:,:,1));
    end
    
    [~, filtered] = findParticles(firstImage,max(max(firstImage)),particleLength,filter,noiseLength);

    f = figure();
    %ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
    set(f, 'Position', [500 300 800 600])
    h1 = subplot(2,2,1);
    set(h1, 'Position', [0.01, 0.5, 0.48, 0.48]);
    imagesc(firstImage)
    ax = gca;
    ax.Visible = 'off';
    h2 = subplot(2,2,2);
    set(h2, 'Position', [0.51, 0.5, 0.48, 0.48]);
    imagesc(filtered)
    ax = gca;
    ax.Visible = 'off';
    
    minThresholdSlider = 0;
    maxThresholdSlider = max(max(filtered));
    
    bgcolor = f.Color;
    
    a = uicontrol('Parent',f,'Style','slider','Units','normalized',...
                'Position',[0.12 0.30 0.77 0.04], ...
                'Tag', 'thresholdSlider', 'UserData', struct('val',0.0, 'image',filtered), ...
                'min',minThresholdSlider, 'max',maxThresholdSlider);
    a.Callback = @(hObject1,ed1) slider_callback1(hObject1,ed1,particleLength,filter,noiseLength);
    al1 = uicontrol('Parent',f,'Style','text','Units','normalized',...
                'Position',[0.08 0.30 0.04 0.03], ...
                'String',sprintf('%.f',minThresholdSlider),'BackgroundColor',bgcolor);
    al2 = uicontrol('Parent',f,'Style','text','Units','normalized',...
                'Position',[0.89 0.30 0.04 0.03], ...
                'String',sprintf('%.f',maxThresholdSlider),'BackgroundColor',bgcolor);
    al3 = uicontrol('Parent',f,'Style','text','Units','normalized',...
                'Position',[0.15 0.36 0.7 0.03], ...
                'String','Threshold','BackgroundColor',bgcolor);
                
    button = uicontrol('Parent', f,'Style','pushbutton','Units','normalized',...
                'Position',[0.35 0.03 0.3 0.1],...
                'Tag', 'okButton', 'UserData', struct('particleLength',particleLength,... 
                'filter', filter, 'noiseLength', noiseLength, 'imageFolder', imageFolder,...
                'uniqueDataNameCounter', uniqueDataNameCounter, 'settings', settings), ...
                'String','Ok');
    button.Callback = @(hObject2,ed2) button_callback(hObject2,ed2, images);
end
function slider_callback1(hObject,eventdata,particleLength,filter,noiseLength)
	sval = hObject.Value;
	hObject.UserData.val = sval;
    adjustThreshold(sval,particleLength,filter,noiseLength)
end
function adjustThreshold(val, particleLength, filter, noiseLength)
    h1 = findobj('Tag','thresholdSlider');
    data = h1.UserData;
    disp('Finding particles... Please wait.')
    [brightPoints, filteredImage] = findParticles(data.image,val,particleLength,filter,noiseLength);
    if(size(brightPoints,1) ~= 0)
        disp([sprintf('%d',size(brightPoints,1)), ' particles found.'])
    else
        disp('No particles found.');
    end
    imagesc(filteredImage)
    hold on
    for i=1:1:size(brightPoints,1)
        viscircles([brightPoints(i,1) brightPoints(i,2)], 1,'Color', 'w', 'LineWidth', 2);
    end
    hold off
    ax = gca;
    ax.Visible = 'off';
end
function button_callback(hObject,eventdata, images)
    h1 = findobj('Tag','thresholdSlider');
    h2 = findobj('Tag','okButton');
    data1 = h1.UserData;
    data2 = h2.UserData;
    %close window
    set(gcf, 'Visible', 'off')
    runMainPTVProgramBlock(data1,data2, images)
end
function runMainPTVProgramBlock(data1,data2, images)
    threshold = data1.val;
    particleLength = data2.particleLength;
    filter = data2.filter;
    noiseLength = data2.noiseLength;
    imageFolder = data2.imageFolder;
    uniqueDataNameCounter = data2.uniqueDataNameCounter;
    settings = data2.settings;
    if(~isempty(settings))
        if(isempty(settings.masks))
            masks = 'No';
        else
            masks = 'Yes';
        end
        if(threshold == 0)
            threshold = num2str(settings.threshold);
        end
        overlap = num2str(settings.overlap);
        subWindow = char(['[',num2str(settings.subWindow),']']);
        searchArea = char(['[',num2str(settings.searchArea),']']);
        trackingSearchRadius = num2str(settings.trackingSearchRadius);
        particleSkippedFrames = num2str(settings.particleSkippedFrames);

        replaceOutliers = settings.replaceOutliers;
        if(replaceOutliers == 2)
            replaceOutliers = char('Without mask');
        elseif(replaceOutliers == 1)
            if(isempty(settings.masks))
                replaceOutliers = char('Without mask');
            else
                replaceOutliers = char('With mask');
            end
        else
            replaceOutliers = char('No');
        end

        ensembleContribution = num2str(settings.ensembleContribution);
        simplifyContribution = num2str(settings.simplifyContribution);
        previousVelocityContribution = num2str(settings.previousVelocityContribution);
        distortedPass = num2str(settings.distortedPass);
        subpixel = char(settings.subpixel);
        trajectoryLength = num2str(settings.trajectoryLength);
        removeLinesOutsideOfMask = char(settings.removeLinesOutsideOfMask);
        
        prompt = {'Use Masks; Yes, No:', 'Particle Brightness Threshold:', 'Particle Diameter in Pixels:', 'Image Filtering; Aggressive, Gaussian, Box:', 'Image Noise Length:', 'PIV Overlap', 'PIV Sub Window: Rectangle of length x and height y; [x,y]', 'PIV Search Area: Rectangle of length x and height y; [x,y]', 'Tracking Search Radius:', 'Number of Frames Particles may Disappear:', 'Replace Outliers; with mask, without mask, no:', 'Ensemble PIV % Contribution: 1 = 100%; Only Ensemble. 0 = 0%; Only Image to Image PIV', 'Simplify Images % Contribution: 1 = 100%; Only Particles Found. 0 = 0%; Using Images As Is.', 'Previous Velocity Field % Contribution: 1 = 100%; Current Velocity Field Is Ignored. 0 = 0%; Only Using Current Velocity Field.', 'Amount of Distorted Passes: Makes the Velocity Field more Accurate; 0,1,2,3', 'Subpixel Accuracy: 3x2, 3x3, 3x3lm, 3x3ls, 5x5lm, 5x5ls or none', 'Plot PIV Signal to Noise Ratio.', 'Remove trajectories smaller than X points:', 'Remove trajectories outside of mask:', 'Plot the Trajectories:'};
        dlg_title = 'PTV inputs';
        num_lines = 1;
        defaultans = {'Yes', num2str(threshold), num2str(particleLength), filter, num2str(noiseLength), overlap, subWindow, searchArea, trackingSearchRadius, particleSkippedFrames, replaceOutliers, ensembleContribution, simplifyContribution, previousVelocityContribution, distortedPass, subpixel, 'Yes', trajectoryLength, removeLinesOutsideOfMask, 'Yes'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    else
        prompt = {'Use Masks; Yes, No:', 'Particle Brightness Threshold:', 'Particle Diameter in Pixels:', 'Image Filtering; Aggressive, Gaussian, Box:', 'Image Noise Length:', 'PIV Overlap', 'PIV Sub Window: Rectangle of length x and height y; [x,y]', 'PIV Search Area: Rectangle of length x and height y; [x,y]', 'Tracking Search Radius:', 'Number of Frames Particles may Disappear:', 'Replace Outliers; with mask, without mask, no:', 'Ensemble PIV % Contribution: 1 = 100%; Only Ensemble. 0 = 0%; Only Image to Image PIV', 'Simplify Images % Contribution: 1 = 100%; Only Particles Found. 0 = 0%; Using Images As Is.', 'Previous Velocity Field % Contribution: 1 = 100%; Current Velocity Field Is Ignored. 0 = 0%; Only Using Current Velocity Field.', 'Amount of Distorted Passes: Makes the Velocity Field more Accurate; 0,1,2,3', 'Subpixel Accuracy: 3x2, 3x3, 3x3lm, 3x3ls, 5x5lm, 5x5ls or none', 'Plot PIV Signal to Noise Ratio.', 'Remove trajectories smaller than X points:', 'Remove trajectories outside of mask:', 'Plot the Trajectories:'};
        dlg_title = 'PTV inputs';
        num_lines = 1;
        defaultans = {'Yes', num2str(threshold), num2str(particleLength), filter, num2str(noiseLength), '0.5', '[30, 30]', '[10, 10]', '2.0', '0', 'Without mask', '0.75', '0.75', '0.5', '3', '3x3ls', 'Yes', '3', 'Yes', 'Yes'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    end
    if(isempty(answer))
        return;
    end
    
    uniqueDataNameCounter = uniqueDataNameCounter + 1;
    assignin('base','uniqueDataNameCounter',uniqueDataNameCounter)
    
    %Resetting to fill up with new settings of this run.
    settings = struct();

    if(strcmp(answer(1), 'y') || strcmp(answer(1), 'Y') || strcmp(answer(1), 'Yes') || strcmp(answer(1), 'yes') || strcmp(answer(1), 'YEs') || strcmp(answer(1), 'YES'))
        if(exist('maskFolder', 'var') && ischar(maskFolder) ~= 0)
            try
                maskFolder = uigetdir(maskFolder, 'Select Mask Folder');
                assignin('base','maskFolder',maskFolder)
            catch
                maskFolder = 0;
                disp('No mask folder selected')
            end
        else
            try
                maskFolder = uigetdir(imageFolder, 'Select Mask Folder');
                assignin('base','maskFolder',maskFolder)
            catch
                maskFolder = 0;
                disp('No mask folder selected')
            end
        end
        if(ischar(maskFolder) ~= 0)
            addpath(maskFolder);
            mask = dir(maskFolder);
            fieldNames = fieldnames(mask);
            masks = rmfield(mask, fieldNames(2:end));
            masks = struct2cell(masks);
            masks(:,1:2) = [];
        else
            masks = [];
        end
    else
        masks = [];
        maskFolder = 0;
    end
    
    settings.images = images;
    settings.imageFolder = imageFolder;
    settings.masks = masks;
    settings.maskFolder = maskFolder;

    threshold = str2double(cell2mat(answer(2)));
    settings.threshold = threshold;
    
    particleLength = str2double(cell2mat(answer(3)));
    settings.particleLength = particleLength;
    
    if(strcmpi(answer(4), 'Aggressive'))
        filter = 'Aggressive';
    elseif(strcmpi(answer(4), 'Gaussian'))
        filter = 'Gaussian';
    elseif(strcmpi(answer(4), 'Box'))
        filter = 'Box';
    else
        filter = 'no';
    end
    settings.filter = filter;
    
    noiseLength = round(str2double(cell2mat(answer(5))));
    settings.noiseLength = noiseLength;
    overlap = str2double(cell2mat(answer(6)));
    settings.overlap = overlap;
    subWindow = regexp(char(answer(7)), '\D*', 'split');
    subWindow = subWindow(~cellfun('isempty',subWindow));
    subWindow = [str2double(subWindow{1}), str2double(subWindow{2})];
    settings.subWindow = subWindow;
    searchArea = regexp(char(answer(8)), '\D*', 'split');
    searchArea = searchArea(~cellfun('isempty',searchArea));
    searchArea = [str2double(searchArea{1}), str2double(searchArea{2})];
    settings.searchArea = searchArea;
    trackingSearchRadius = str2double(cell2mat(answer(9)));
    settings.trackingSearchRadius = trackingSearchRadius;
    particleSkippedFrames = str2double(cell2mat(answer(10)));
    settings.particleSkippedFrames = particleSkippedFrames;
    
    if(strcmpi(answer(11), 'Without mask'))
        replaceOutliers = 2;
    elseif(strcmpi(answer(11), 'With mask'))
        if(isempty(masks))
            replaceOutliers = 2;
        else
            replaceOutliers = 1;
        end
    else
        replaceOutliers = 0;
    end
    settings.replaceOutliers = replaceOutliers;
    
    ensembleContribution = str2double(cell2mat(answer(12)));
    settings.ensembleContribution = ensembleContribution;
    
    simplifyContribution = str2double(cell2mat(answer(13)));
    settings.simplifyContribution = simplifyContribution;
    if(simplifyContribution > 0 && simplifyContribution <= 1)
        prompt = {'Blur image (improved PIV); none, gaussian:', 'Blur Strength:'};
        dlg_title = 'Simplify Image inputs. This creates white dots where the particles have been found to help with the PIV. A little blur makes the dots more realistic as particles.';
        num_lines = 1;
        defaultans = {'gaussian', '0.5'};
        drawAnswer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if(strcmpi(drawAnswer(1), 'none'))
            blur = 'none';
        else
            blur = 'gaussian';
        end
        sigma = str2double(cell2mat(drawAnswer(2)));
    elseif(simplifyContribution == 0)
        blur = 'none';
        sigma = 0.5;
    end
    settings.blur = blur;
    settings.sigma = sigma;
    
    previousVelocityContribution = str2double(cell2mat(answer(14)));
    settings.previousVelocityContribution = previousVelocityContribution;
    
    distortedPass = str2double(cell2mat(answer(15)));
    settings.distortedPass = distortedPass;
    
    subpixel = answer(16);
    settings.subpixel = subpixel;
    
    if(strcmpi(subpixel,'3x2'))
        subpixel = @subpixel3x2;
    elseif(strcmpi(subpixel,'3x3'))
        subpixel = @subpixel3x3;
    elseif(strcmpi(subpixel,'3x3lm'))
        subpixel = @subpixel3x3lm;
    elseif(strcmpi(subpixel,'3x3ls'))
        subpixel = @subpixel3x3ls;
    elseif(strcmpi(subpixel,'5x5lm'))
        subpixel = @subpixel5x5lm;
    elseif(strcmpi(subpixel,'5x5ls'))
        subpixel = @subpixel5x5ls;
    else
        subpixel = @subpixelnone;
    end
    
    trajectoryLength = str2double(cell2mat(answer(18)));
    settings.trajectoryLength = trajectoryLength;
    removeLinesOutsideOfMask = answer(19);
    settings.removeLinesOutsideOfMask = removeLinesOutsideOfMask;
    
    assignin('base',sprintf('settings_%d',uniqueDataNameCounter),settings)
    
    [firstFoundTrackingPoints, piv, Uint, Vint, ensembleUint, ensembleVint] = PTV(images, masks, threshold, particleLength, filter, noiseLength, overlap, subWindow, searchArea, distortedPass, trackingSearchRadius, particleSkippedFrames, replaceOutliers, ensembleContribution, simplifyContribution, previousVelocityContribution, blur, sigma, subpixel);
    assignin('base',sprintf('firstFoundTrackingPoints_%d',uniqueDataNameCounter),firstFoundTrackingPoints)
    assignin('base',sprintf('piv_%d',uniqueDataNameCounter),piv)
    assignin('base',sprintf('Uint_%d',uniqueDataNameCounter),Uint)
    assignin('base',sprintf('Vint_%d',uniqueDataNameCounter),Vint)
    assignin('base',sprintf('ensembleUint_%d',uniqueDataNameCounter),ensembleUint)
    assignin('base',sprintf('ensembleVint_%d',uniqueDataNameCounter),ensembleVint)
    
    if(strcmpi(answer(17), 'y') || strcmpi(answer(17), 'Ye') || strcmpi(answer(17), 'Yes'))
        fontSize = 30;
        figure
        imagesc(piv.snr);
        titleText = 'PIV Signal to Noise Ratio.';
        title(titleText,'FontSize',fontSize)
        set(gca, 'FontSize', fontSize)
    end
    
    %Removing short tracks.
    if(trajectoryLength > 0)
        sortedFFTP = sortrows(firstFoundTrackingPoints, 4);
        findDifferenceIndex = find(diff(sortedFFTP(:,4)) ~= 0);
        if(~isempty(findDifferenceIndex))
            trajectoryLengthsInArray = diff(findDifferenceIndex);
            if(findDifferenceIndex(1) <= trajectoryLength)
                sortedFFTP(1:findDifferenceIndex(1),:) = 0;
            end
            if(size(sortedFFTP,1) - findDifferenceIndex(end) <= trajectoryLength)
                sortedFFTP(findDifferenceIndex(end)+1:size(sortedFFTP,1),:) = 0;
            end
            for i=1:size(trajectoryLengthsInArray,1)-1
                if(trajectoryLengthsInArray(i) <= trajectoryLength)
                    sortedFFTP(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),:) = 0;
                end
            end

            %Removing rows of zeros.
            sortedFFTP( ~any(sortedFFTP,2), : ) = [];
        end
    else
        sortedFFTP = sortrows(firstFoundTrackingPoints, 4);
    end
    
    if(strcmpi(answer(19), 'y') || strcmpi(answer(19), 'Ye') || strcmpi(answer(19), 'Yes'))
        fontSize = 30;
        if(isempty(masks))
            mask = 0;
        else
            mask = imread(masks{1});
        end
        %Removing lines that cross over the mask. The coordinates may be within
        %the mask, but the line between them is outside of the mask.
        sortedFFTP = RemoveLinesOutsideOfMask(sortedFFTP, mask, 1000);
    end
    
    assignin('base',sprintf('sortedFFTP_%d',uniqueDataNameCounter),sortedFFTP)
    
    if(strcmpi(answer(20), 'y') || strcmpi(answer(20), 'Ye') || strcmpi(answer(20), 'Yes'))
        if(isempty(masks))
            backgroundMask = 0;
        else
            mask = imread(masks{1});
            backgroundMask = bwdist(imcomplement(mask));
        end
        if(ensembleContribution ~= 0)
            backgroundVelocity = sqrt(ensembleUint.^2+ensembleVint.^2);
        else
            backgroundVelocity = sqrt(Uint.^2+Vint.^2);
        end
        image = imread(images{1});
        if(size(image,3) == 3)
            backgroundImage = rgb2gray(image);
        else
            backgroundImage = image;
        end
        
        %This makes it so that we skip plotting for mask as background if
        %there are no selected masks.
        if(isempty(masks))
            startNr = 2;
        else
            startNr = 1;
        end
        for selectedBackground=startNr:3
            figure
            markerSize = round(fontSize*2.5);
            if(selectedBackground==1)
            	contourf(backgroundMask)
                c = colorbar('southoutside','FontSize',fontSize);
                c.Label.String = 'Distance to Vein Margin in Pixels';
                titleText = ['Particle to Border Distance. Threshold: ', num2str(threshold), '.'];
            elseif(selectedBackground==2)
                contourf(backgroundVelocity)
                c = colorbar('southoutside','FontSize',fontSize);
                c.Label.String = 'Velocity in Pixels per Frame';
                titleText = ['Particle Trajectories on Magnitude of Velocity Field. Threshold: ', num2str(threshold), '.'];
            elseif(selectedBackground==3)
                imagesc(backgroundImage)
                titleText = ['Particle Trajectories on Source Image. Threshold: ', num2str(threshold), '.'];
            end
            
            set(gca, 'YDir', 'reverse')
            set(gca, 'FontSize', fontSize);
            title(titleText,'FontSize',fontSize)
            xlabel('Pixels','FontSize',fontSize) % x-axis label
            ylabel('Pixels','FontSize',fontSize) % y-axis label
            hold on

            for i=1:size(sortedFFTP,1)-1
                if(sortedFFTP(i,4) == sortedFFTP(i+1,4))
                    hold on
                    line([sortedFFTP(i,1) sortedFFTP(i+1,1)], [sortedFFTP(i,2) sortedFFTP(i+1,2)], 'Color', 'w', 'linewidth', 2.5)
                    if(i ~= size(sortedFFTP,1)-1)
                        if(sortedFFTP(i+1,4) ~= sortedFFTP(i+2,4))
                            hold on
                            scatter(sortedFFTP(i+1,1), sortedFFTP(i+1,2),markerSize, 'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth', 1.5)
                        end
                    else
                        hold on
                        scatter(sortedFFTP(i+1,1), sortedFFTP(i+1,2),markerSize, 'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth', 1.5)
                    end
                end
            end
        end
    end
end