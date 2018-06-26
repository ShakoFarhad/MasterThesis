function []=particleAverageVelocityQuiverPlot(FFTP, bgImg)
    trajectoryEndPoints = sortrows(FFTP,4);
    findDifferenceIndexID = find(diff(trajectoryEndPoints(:,4)) ~= 0);
    trajectoryLengthsInArray = diff(findDifferenceIndexID);
    if(findDifferenceIndexID(1) >= 2)
        trajectoryEndPoints(1:findDifferenceIndexID(1)-1,:) = 0;
    end
    if(size(trajectoryEndPoints,1) - findDifferenceIndexID(end) >= 2)
        trajectoryEndPoints(findDifferenceIndexID(end)+1:size(trajectoryEndPoints,1)-1,:) = 0;
    end
    for i=1:size(trajectoryLengthsInArray,1)
        if(trajectoryLengthsInArray(i) >= 2)
            trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1)-1,:) = 0;
        else
            trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1),:) = 0;
        end
    end
    trajectoryEndPoints(trajectoryEndPoints(:,4)==0,:) = [];

    averageTV = averageTrackVelocity(FFTP);

    fontSize = 30;
    figure
    colormap(winter)
    contourf(bgImg)
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', fontSize);
    hold on;
    quiver(trajectoryEndPoints(:,1),trajectoryEndPoints(:,2),averageTV(:,1),averageTV(:,2), 0, 'w', 'linewidth', 1.6)
    hold off
    titleText = 'Magnitude of Velocity Field';
    title(titleText,'FontSize',fontSize)
    c = colorbar('southoutside','FontSize',fontSize);
    c.Label.String = 'Velocity in Pixels per Frame';
    xlabel('Pixels','FontSize',fontSize) % x-axis label
    ylabel('Pixels','FontSize',fontSize) % y-axis label
end