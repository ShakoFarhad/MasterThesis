prompt = {'Outer Trajectories: 0, 1, 2 etc; 0 = no data.', 'Inner Trajectories: 0, 1, 2 etc; 0 = no data.', 'Full Trajectories: 0, 1, 2 etc; 0 = no data.'};
dlg_title = 'Plot Trajectories.';
num_lines = 1;
if(exist('uniqueDataNameCounter', 'var'))
    if(uniqueDataNameCounter == 1)
        defaultans = {num2str(uniqueDataNameCounter), '0', '0'};
    elseif(uniqueDataNameCounter == 2)
        defaultans = {num2str(uniqueDataNameCounter-1),num2str(uniqueDataNameCounter), '0'};
    elseif(uniqueDataNameCounter == 3)
        defaultans = {num2str(uniqueDataNameCounter-2),num2str(uniqueDataNameCounter-1), num2str(uniqueDataNameCounter)};
    else
        defaultans = {num2str(uniqueDataNameCounter-2),num2str(uniqueDataNameCounter-1), num2str(uniqueDataNameCounter)};
    end
else
    defaultans = {'0','0','0'};
end
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if(isempty(answer))
    outerData = 0;
    innerData = 0;
    fullData = 0;
else
    outerData = str2double(cell2mat(answer(1)));
    innerData = str2double(cell2mat(answer(2)));
    fullData = str2double(cell2mat(answer(3)));
end

%Input values used for calculation of the particle trajectories's angle to
%a reference line
framesMovement = 1;
useMaskLine = 1;
angleFull = 10;
angleCombined = 10;

fontSize = 30;

if(outerData ~= 0)
    %Outer Data
    try
        firstFoundTrackingPointsOuter = eval(sprintf('firstFoundTrackingPoints_%d',outerData));
    catch
        disp([sprintf('firstFoundTrackingPoints_%d',outerData) ' does not exist. Please run runPTV script.'])
        return
    end
    settingsOuter = eval(sprintf('settings_%d',outerData));
    sortedFFTPOuter = eval(sprintf('sortedFFTP_%d',outerData));
    UintOuter = eval(sprintf('Uint_%d',outerData));
    VintOuter = eval(sprintf('Vint_%d',outerData));
    ensembleUintOuter = eval(sprintf('ensembleUint_%d',outerData));
    ensembleVintOuter = eval(sprintf('ensembleVint_%d',outerData));
end

if(innerData ~= 0)
    %Inner Data
    try
        firstFoundTrackingPointsInner = eval(sprintf('firstFoundTrackingPoints_%d',innerData));
    catch
        disp([sprintf('firstFoundTrackingPoints_%d',innerData) ' does not exist. Please run runPTV script.'])
        return
    end
    settingsInner = eval(sprintf('settings_%d',innerData));
    sortedFFTPInner = eval(sprintf('sortedFFTP_%d',innerData));
    UintInner = eval(sprintf('Uint_%d',innerData));
    VintInner = eval(sprintf('Vint_%d',innerData));
    ensembleUintInner = eval(sprintf('ensembleUint_%d',innerData));
    ensembleVintInner = eval(sprintf('ensembleVint_%d',innerData));
end

if(fullData ~= 0)
    %Full Data
    try
        firstFoundTrackingPointsFull = eval(sprintf('firstFoundTrackingPoints_%d',fullData));
    catch
        disp([sprintf('firstFoundTrackingPoints_%d',fullData) ' does not exist. Please run runPTV script.'])
        return
    end
    settingsFull = eval(sprintf('settings_%d',fullData));
    sortedFFTPFull = eval(sprintf('sortedFFTP_%d',fullData));
    UintFull = eval(sprintf('Uint_%d',fullData));
    VintFull = eval(sprintf('Vint_%d',fullData));
    ensembleUintFull = eval(sprintf('ensembleUint_%d',fullData));
    ensembleVintFull = eval(sprintf('ensembleVint_%d',fullData));
    magnitudeVelocityFull = sqrt(UintFull.^2+VintFull.^2);

    masksFull = settingsFull.masks;
    maskFolderFull = settingsFull.maskFolder;
    maskFull = imread(char(strcat(maskFolderFull,'\', masksFull)));
    distanceImageFull = bwdist(imcomplement(maskFull));

    dataTypeFull = 'Full';
    fitPowerLaw(firstFoundTrackingPointsFull, distanceImageFull, dataTypeFull)

    [thetaFull, signChangeCountFull, trajectoryEndPointsFull, totalAmountOfTrajectoriesFull] = CalculateParticleAngle(sortedFFTPFull, maskFull, framesMovement, useMaskLine);
    
    angleFull
    framesMovement
    signChangeCountFull
    percentageOfTrajectoriesThetaFull = size(thetaFull,1)/totalAmountOfTrajectoriesFull*100
    percentageOfTrajectoriesAngleFull = sum(thetaFull>=angleFull)/totalAmountOfTrajectoriesFull*100
    thetaFullMinMax = [min(thetaFull), max(thetaFull)]
    
    figure
    contourf(magnitudeVelocityFull)    
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', fontSize);
    titleText = 'Magnitude of Velocity Field';
    title(titleText,'FontSize',fontSize)
    c = colorbar('southoutside','FontSize',fontSize);
    c.Label.String = 'Velocity in Pixels per Frame';
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label

    particleAverageVelocityQuiverPlot(sortedFFTPFull, magnitudeVelocityFull)
    particleAverageVelocityQuiverPlot(firstFoundTrackingPointsFull, magnitudeVelocityFull)
    
    
    threshold = settingsFull.threshold;
    particleLength = settingsFull.particleLength;
    filter = settingsFull.filter;
    noiseLength = settingsFull.noiseLength;
    images = settingsFull.images;
    
    plotNumberOfParticlesFoundPerImage(images, threshold, particleLength, filter, noiseLength)
end

if(innerData ~= 0 && outerData ~= 0 && fullData ~= 0)
    dataTypeInner = 'Inner';
    fitPowerLaw(firstFoundTrackingPointsInner, distanceImageFull, dataTypeInner)
    
    %Combining Outer and Inner
    UintCombined = UintOuter + UintInner;
    VintCombined = VintOuter + VintInner;
    ensembleUintCombined = ensembleUintOuter + ensembleUintInner;
    ensembleVintCombined = ensembleVintOuter + ensembleVintInner;
    magnitudeVelocityFieldCombined = sqrt(UintCombined.^2 + VintCombined.^2);

    thresholdOuter = settingsOuter.threshold;
    thresholdInner = settingsInner.threshold;
    
    dataTypeCombined = 'Outer and Inner';
    fitPowerLaw(firstFoundTrackingPointsCombined, distanceImageFull, dataTypeCombined)

    [thetaOuter, signChangeCountOuter, trajectoryEndPointsOuter, totalAmountOfTrajectoriesOuter] = CalculateParticleAngle(sortedFFTPOuter, maskFull, framesMovement, useMaskLine);
    [thetaInner, signChangeCountInner, trajectoryEndPointsInner, totalAmountOfTrajectoriesInner] = CalculateParticleAngle(sortedFFTPInner, maskFull, framesMovement, useMaskLine);
    
    angleCombined = angleCombined
    framesMovement
    percentageOfTrajectoriesAngleOuter = sum(thetaOuter>=angleCombined)/totalAmountOfTrajectoriesOuter*100
    percentageOfTrajectoriesAngleInner = sum(thetaInner>=angleCombined)/size(trajectoryEndPointsInner,1)*100
    percentageOfTrajectoriesAngleCombined = (sum(thetaOuter>=angleCombined) + sum(thetaInner>=angleCombined))/(totalAmountOfTrajectoriesInner + totalAmountOfTrajectoriesOuter)*100
    
    percentageOfTrajectoriesThetaOuter = size(thetaOuter,1)/totalAmountOfTrajectoriesOuter*100
    percentageOfTrajectoriesThetaInner = size(thetaInner,1)/totalAmountOfTrajectoriesInner*100
    percentageOfTrajectoriesThetaCombined = (size(thetaInner,1) + size(thetaOuter,1))/(totalAmountOfTrajectoriesInner + totalAmountOfTrajectoriesOuter)*100
    
    thetaOuterMinMax = [min(thetaOuter), max(thetaOuter)]
    thetaInnerMinMax = [min(thetaInner), max(thetaInner)]
    
    averageTPOuter = averageTrackPosition(sortedFFTPOuter);
    averageTVOuter = averageTrackVelocity(sortedFFTPOuter);
    averageTPInner = averageTrackPosition(sortedFFTPInner);
    averageTVInner = averageTrackVelocity(sortedFFTPInner);

    figure
    contourf(magnitudeVelocityFieldCombined)    
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', fontSize);
    titleText = 'Magnitude of Velocity Field';
    title(titleText,'FontSize',fontSize)
    c = colorbar('southoutside','FontSize',fontSize);
    c.Label.String = 'Velocity in Pixels per Frame';
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label
    
    figure
    contourf(magnitudeVelocityFieldCombined)
    colormap(winter)
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', fontSize);
    hold on;
    quiver(averageTPOuter(:,1),averageTPOuter(:,2),averageTVOuter(:,1),averageTVOuter(:,2), 0, 'y')
    hold on
    quiver(averageTPInner(:,1),averageTPInner(:,2),averageTVInner(:,1),averageTVInner(:,2), 0, 'y')
    titleText = 'Magnitude of Velocity Field';
    title(titleText,'FontSize',fontSize)
    c = colorbar('southoutside','FontSize',fontSize);
    c.Label.String = 'Velocity in Pixels per Frame';
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label

    averageTPOuter = averageTrackPosition(firstFoundTrackingPointsOuter);
    averageTVOuter = averageTrackVelocity(firstFoundTrackingPointsOuter);
    averageTPInner = averageTrackPosition(firstFoundTrackingPointsInner);
    averageTVInner = averageTrackVelocity(firstFoundTrackingPointsInner);
    
    figure
    contourf(magnitudeVelocityFieldCombined)
    colormap(winter)
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', fontSize);
    hold on;
    quiver(averageTPOuter(:,1),averageTPOuter(:,2),averageTVOuter(:,1),averageTVOuter(:,2), 0, 'y')
    hold on
    quiver(averageTPInner(:,1),averageTPInner(:,2),averageTVInner(:,1),averageTVInner(:,2), 0, 'y')
    titleText = 'Magnitude of Velocity Field';
    title(titleText,'FontSize',fontSize)
    c = colorbar('southoutside','FontSize',fontSize);
    c.Label.String = 'Velocity in Pixels per Frame';
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label

    averageTPOuter = round(averageTrackPosition(firstFoundTrackingPointsOuter));
    im1 = imread(settingsFull.images{1});
    if(size(im1,3) == 3)
        im1 = im2double(rgb2gray(im1));
    else
        im1 = im2double(im1(:,:,1));
    end
    particleLength = settingsFull.particleLength;
    blur = settingsFull.blur;
    sigma = settingsFull.sigma;

    ensembleContribution = settingsOuter.ensembleContribution;

    if(ensembleContribution ~= 0)
        backgroundVelocity = sqrt(ensembleUintCombined.^2+ensembleVintCombined.^2);
    else
        backgroundVelocity = sqrt(UintCombined.^2+VintCombined.^2);
    end
    image = imread(images{1});
    if(size(image,3) == 3)
        backgroundImage = rgb2gray(image);
    else
        backgroundImage = image;
    end

    for selectedBackground=1:3
        figure
        fontSize = 30;
        markerSize = round(fontSize*2.5);
        if(selectedBackground==1)
            contourf(distanceImageFull)
            c = colorbar('southoutside','FontSize',fontSize);
            c.Label.String = 'Distance to Vein Margin in Pixels';
            titleText = ['Particle to Border Distance. Threshold Outer and Inner: ', num2str(thresholdOuter), ', ', num2str(thresholdInner), '.'];
        elseif(selectedBackground==2)
            contourf(backgroundVelocity)
            c = colorbar('southoutside','FontSize',fontSize);
            c.Label.String = 'Velocity in Pixels per Frame';
            titleText = ['Particle Trajectories on Velocity Field. Threshold Outer and Inner: ', num2str(thresholdOuter), ', ', num2str(thresholdInner), '.'];
        elseif(selectedBackground==3)
            contourf(backgroundImage)
            titleText = ['Particle Trajectories on Source Image. Threshold Outer and Inner: ', num2str(thresholdOuter), ', ', num2str(thresholdInner), '.'];
        end

        set(gca, 'YDir', 'reverse')
        set(gca, 'FontSize', fontSize);
        title(titleText,'FontSize',fontSize)
        xlabel('Pixels','FontSize',fontSize) % x-axis label
        ylabel('Pixels','FontSize',fontSize) % y-axis label
        hold on

        for i=1:size(sortedFFTPOuter,1)-1
            if(sortedFFTPOuter(i,4) == sortedFFTPOuter(i+1,4))
                hold on
                line([sortedFFTPOuter(i,1) sortedFFTPOuter(i+1,1)], [sortedFFTPOuter(i,2) sortedFFTPOuter(i+1,2)], 'Color', 'w', 'linewidth', 2.5)
                if(i ~= size(sortedFFTPOuter,1)-1)
                    if(sortedFFTPOuter(i+1,4) ~= sortedFFTPOuter(i+2,4))
                        hold on
                        scatter(sortedFFTPOuter(i+1,1), sortedFFTPOuter(i+1,2),markerSize, 'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth', 1.5)
                    end
                end
            end
        end
        for i=1:size(sortedFFTPInner,1)-1
            if(sortedFFTPInner(i,4) == sortedFFTPInner(i+1,4))
                hold on
                line([sortedFFTPInner(i,1) sortedFFTPInner(i+1,1)], [sortedFFTPInner(i,2) sortedFFTPInner(i+1,2)], 'Color', 'w', 'linewidth', 2.5)
                if(i ~= size(sortedFFTPInner,1)-1)
                    if(sortedFFTPInner(i+1,4) ~= sortedFFTPInner(i+2,4))
                        hold on
                        scatter(sortedFFTPInner(i+1,1), sortedFFTPInner(i+1,2),markerSize, 'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth', 1.5)
                    end
                end
            end
        end
    end
end