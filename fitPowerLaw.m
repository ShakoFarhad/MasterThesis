function []=fitPowerLaw(firstFoundTrackingPoints, distanceImage, dataType, binWidth, binStart, binEnd)
    % Calculate Particle Angle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % firstFoundTrackingPoints: Comes from RunPTV.m
    % distanceImage: Distances from every pixel to the border of the mask
    % histBinStart: Histogram Bin Start
    % histBinEnd: Histogram Bin End
    % histBinWidth: Histogram Bin Width
    narginchk(2,5)

    particleDistances = [diag(distanceImage(firstFoundTrackingPoints(:,2), firstFoundTrackingPoints(:,1))), firstFoundTrackingPoints(:,4)];

    IDsortedTracks = sortrows(particleDistances,2);
    findDifferenceIndex = find(diff(IDsortedTracks(:,2)) ~= 0);
    trajectoryLengthsInArray = diff(findDifferenceIndex);
    
    meanParticleDistances = zeros(size(trajectoryLengthsInArray,1)+2,1);
    meanD = mean(IDsortedTracks(1:findDifferenceIndex(1),1));
    meanParticleDistances(1,:) = meanD;
    for i=1:size(trajectoryLengthsInArray,1)
        meanD = mean(IDsortedTracks(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),1));
        meanParticleDistances(i+1,:) = meanD;
    end
    meanD = mean(IDsortedTracks(findDifferenceIndex(end)+1:size(IDsortedTracks,1),1));
    meanParticleDistances(end,:) = meanD;
    
    if ~exist('binWidth', 'var')
        binWidth = 1;
    end
    if ~exist('binStart', 'var')
        binStart = 0;
    end
    if ~exist('binEnd', 'var')
        binEnd = double(max(meanParticleDistances));
    end

    fontSize = 30;
    figure
    binCtrs = binStart:binWidth:binEnd;
    hist(meanParticleDistances(:,1), binCtrs)
    h_gca = gca;
    h = h_gca.Children;
    h.FaceColor = [.91 .91 .91];
    h.EdgeColor = [.98 .98 .98];
    counts = hist(meanParticleDistances(:,1),binCtrs);

    t = binCtrs;
    y = counts;
    F = @(x)(x(2))*((t).^(x(1))).*((binEnd-t).^(x(1)))-y;
    x0 = [binStart+1 binEnd-1];
    x = lsqnonlin(F,x0);

    hold on
    plot(t,y,'bo', t, F(x)+y, 'r-', 'LineWidth', 1.5);
    legend('Histogram of Data', 'Histogram Centers', 'Fitted Power Law to Data')
    titleText = [sprintf('Fitted b x^a (%.f-x)^a, a = %.4f', binEnd, x(1)), sprintf(', b = %.4f', x(2)), '. ', dataType, ' Particle Trajectories'];
    title(titleText)
    axis([binStart-1 binEnd+1 -inf 1.03*max(y)])
    xlabel('Distance from particle to boundary in pixels') % x-axis label
    ylabel('Number of particle trajectories') % y-axis label
    set(gca, 'FontSize', fontSize)
end