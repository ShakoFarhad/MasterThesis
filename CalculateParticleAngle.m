function [varargout]=CalculateParticleAngle(FFTP, mask, framesMovement, useMaskLine)
    % Calculate Particle Angle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sortedFFTP: Comes from RunPTV.m
    % mask: A binary mask image
    % angle: Angle between lines in degrees
    % framesMovement: How many frames the particle moves in the same 
    % direction with the same velocity
    % useMaskLine: Draw lines on vein or find lines through mask; 1 or 0

    %Minimum 2 inputs and maximum 4 inputs
    narginchk(2,4)
    if ~exist('framesMovement', 'var')
        framesMovement = 1;
    end
    
    if ~exist('useMaskLine', 'var')
        useMaskLine = 1;
    end

    fontSize = 30;
    if(useMaskLine)
        if(size(mask,3) == 3)
            maskDistances = bwdist(rgb2gray(imcomplement(mask)));
        else
            maskDistances = bwdist(imcomplement(mask));
        end
        [maxDistances, y] = max(maskDistances);
        x = 1:size(y,2);
    else
        figure
        imshow(im2double(mask(:,:,1)))
        title('Trace the vein by clicking on the image from left to right. Press Enter to exit.','FontSize',fontSize)
        xlabel('Pixels','FontSize',fontSize) % x-axis label
        ylabel('Pixels','FontSize',fontSize) % y-axis label

        %Amount of points to click on the picture to define the line.
        [x,y] = ginput;
        x = transpose(x);
        y = transpose(y);
        close all
    end

    figure
    imshow(mask)
    hold on
    plot(x,y,'b-', 'LineWidth', 5)
    hold on
    scatterPlotMarkerSize = floor(0.3*(size(mask,1) + size(mask,2))/2);
    scatter(x,y, scatterPlotMarkerSize, 'MarkerFaceColor','k','MarkerEdgeColor','r', 'LineWidth', 3)
    title('The selected endpoints.','FontSize',fontSize)
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label

    framesMovement = 1/framesMovement;
    %Finding the endpoints for the different trajectories
    %If the velocity is greater than its velocity times framesMovement then 
    %we will keep that track because then it is likely that it would keep
    %moving into the vein wall.
    totalAmountOfTrajectories = 0;
    trajectoryEndPoints = sortrows(FFTP,4);
    findDifferenceIndexID = find(diff(trajectoryEndPoints(:,4)) ~= 0);
    trajectoryLengthsInArray = diff(findDifferenceIndexID);
    if(findDifferenceIndexID(1) >= 2)
        distanceFromVeinWall = maskDistances(trajectoryEndPoints(findDifferenceIndexID(1),2),trajectoryEndPoints(findDifferenceIndexID(1),1));
        if(trajectoryEndPoints(findDifferenceIndexID(1),6)^2 + trajectoryEndPoints(findDifferenceIndexID(1),7)^2 >= (framesMovement*distanceFromVeinWall)^2)
            trajectoryEndPoints(1:findDifferenceIndexID(1)-2,:) = 0;
            totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
        else
            trajectoryEndPoints(1:findDifferenceIndexID(1),:) = 0;
            totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
        end
    else
        trajectoryEndPoints(1:findDifferenceIndexID(1),:) = 0;
    end
    if(size(trajectoryEndPoints,1) - findDifferenceIndexID(end) >= 2)
        distanceFromVeinWall = maskDistances(trajectoryEndPoints(size(trajectoryEndPoints,1),2),trajectoryEndPoints(size(trajectoryEndPoints,1),1));
        if(trajectoryEndPoints(size(trajectoryEndPoints,1),6)^2 + trajectoryEndPoints(size(trajectoryEndPoints,1),7)^2 >= (framesMovement*distanceFromVeinWall)^2)
            trajectoryEndPoints(findDifferenceIndexID(end)+1:size(trajectoryEndPoints,1)-2,:) = 0;
            totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
        else
            trajectoryEndPoints(findDifferenceIndexID(end)+1:size(trajectoryEndPoints,1),:) = 0;
            totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
        end
    else
        trajectoryEndPoints(findDifferenceIndexID(end)+1:size(trajectoryEndPoints,1),:) = 0;
    end
    for i=1:size(trajectoryLengthsInArray,1)
        distanceFromVeinWall = maskDistances(trajectoryEndPoints(findDifferenceIndexID(i+1),2),trajectoryEndPoints(findDifferenceIndexID(i+1),1));
        if(trajectoryLengthsInArray(i) >= 2)
            if(trajectoryEndPoints(findDifferenceIndexID(i+1),6)^2 + trajectoryEndPoints(findDifferenceIndexID(i+1),7)^2 >= (framesMovement*distanceFromVeinWall)^2)
                trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1)-2,:) = 0;
                totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
            else
                trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1),:) = 0;
                totalAmountOfTrajectories = totalAmountOfTrajectories + 1;
            end
        else
            trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1),:) = 0;
        end
    end
    trajectoryEndPoints(trajectoryEndPoints(:,4)==0,:) = [];

    %Calculating the slopes from the endpoints.
    theta = zeros(size(trajectoryEndPoints,1),1);
    for j=1:size(x,2)-1
        if(j==1)
            if(x(j) < 1)
                x(j) = 1;
            end
            if(y(j) < 1)
                y(j) = 1;
            end
            borderX1 = 1;
            borderX2 = x(j+1);
        elseif(j==size(x,1)-1)
            if(x(j+1) > size(mask,2))
                x(j+1) = size(mask,2);
            end
            if(y(j+1) > size(mask,1))
                y(j+1) = size(mask,1);
            end
            borderX1 = x(j);
            borderX2 = size(mask,2);
        else
            if(x(j) < 1)
                x(j) = 1;
            end
            if(x(j+1) > size(mask,2))
                x(j+1) = size(mask,2);
            end
            if(y(j) < 1)
                y(j) = 1;
            end
            if(y(j+1) > size(mask,1))
                y(j+1) = size(mask,1);
            end
            borderX1 = x(j);
            borderX2 = x(j+1);
        end
        %Calculating the slope for the line that was drawn.
        mainLineSlope = (y(j+1)-y(j))/(x(j+1)-x(j));
        theta1 = radtodeg(atan(mainLineSlope));
        for i=1:2:size(trajectoryEndPoints,1)-1
            if((trajectoryEndPoints(i,1) >= borderX1) && (trajectoryEndPoints(i,1) < borderX2))
                if(y(j) == y(j+1))
                    yLine = y(j);
                else
                    mainLineBias = y(j) - mainLineSlope*x(j);
                    yLine = mainLineSlope*trajectoryEndPoints(i,1) + mainLineBias;
                end
                %Checking if particle is above or below the drawn line/main
                %line.
                if(trajectoryEndPoints(i,2) >= yLine)
                    %Checking if the particle is moving up towards the wall
                    if(trajectoryEndPoints(i,7) <= 0)
                        if(trajectoryEndPoints(i+1,1) - trajectoryEndPoints(i,1) == 0)
                            theta(i) = 90-theta1;
                        else
                            %Calculating the angle between the drawn line and the trajectory lines.
                            trajectorySlope = (trajectoryEndPoints(i+1,2) - trajectoryEndPoints(i,2))/(trajectoryEndPoints(i+1,1) - trajectoryEndPoints(i,1));
                            theta(i) = radtodeg(atan(abs((mainLineSlope - trajectorySlope)/(1+trajectorySlope*mainLineSlope))));
                        end
                    end
                else
                    %Checking if the particle is moving down towards the wall
                    if(trajectoryEndPoints(i,7) >= 0)
                        if(trajectoryEndPoints(i+1,1) - trajectoryEndPoints(i,1) == 0)
                            theta(i) = 90-theta1;
                        else
                            %Calculating the angle between the drawn line and the trajectory lines.
                            trajectorySlope = (trajectoryEndPoints(i+1,2) - trajectoryEndPoints(i,2))/(trajectoryEndPoints(i+1,1) - trajectoryEndPoints(i,1));
                            theta(i) = radtodeg(atan(abs((mainLineSlope - trajectorySlope)/(1+trajectorySlope*mainLineSlope))));
                        end
                    end
                end
            end
        end
    end
    theta(theta(:)==0) = [];
    
    %Calculating the complexity of the vein; how many local maxima and
    %minima it has. First we derivate, then we remove the zeros as these
    %are plateaus, and then we count every sign change.
    pointsBeforeSignChangeCount = 0;
    amountOfPointsBeforeSignChange = 3;
    numberOfPixelsDifferenceAtSignChangeLowerBound = 2;
    numberOfPixelsDifferenceAtSignChangeUpperBound = round(size(mask,1)/3);
    signChangeCount = 0;
    yDiff = diff(y);
    yDiff(yDiff == 0) = [];
    
    for i=1:size(yDiff,2)-1
        pointsBeforeSignChangeCount = pointsBeforeSignChangeCount + 1;
        if(sign(yDiff(i)) ~= sign(yDiff(i+1)))
            %This check is to make sure that we do not count small 1 pixel
            %local minima or maxima. These may be just noise.
            if((abs(yDiff(i+1)) >= numberOfPixelsDifferenceAtSignChangeLowerBound && abs(yDiff(i+1)) >= numberOfPixelsDifferenceAtSignChangeUpperBound) || pointsBeforeSignChangeCount >= amountOfPointsBeforeSignChange)
                signChangeCount = signChangeCount + 1;
            end
            pointsBeforeSignChangeCount = 0;
        end
    end
    
    
    %return
    if nargout == 1
        varargout{1} = theta;
    elseif nargout == 2
        varargout{1} = theta;
        varargout{2} = signChangeCount;
    elseif nargout == 3
        varargout{1} = theta;
        varargout{2} = signChangeCount;
        varargout{3} = trajectoryEndPoints;
    elseif nargout == 4
        varargout{1} = theta;
        varargout{2} = signChangeCount;
        varargout{3} = trajectoryEndPoints;
        varargout{4} = totalAmountOfTrajectories;
    end
end