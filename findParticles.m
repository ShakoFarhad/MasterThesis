function coordinates=findParticles(image,threshold,particleLength,filter, noiseLength, visualize)
%Finds local maxima in an image to pixel level accuracy.
%
%coordinates=findParticles(image,threshold,particleLength,filter,noiseLength)
%
%image - Image of particles.
%threshold - The brightness cutoff point.
%particleLength - The diameter of a particle in pixels.
%filter - Three options. 'bandpass', 'lowpass' and 'no' for no filter.
%noiseLength - The diameter of noise particles in pixels.
%visualize - 'plot' plots dots on the original image.
    %Minimum 1 inputs and maximum 6 inputs
    narginchk(1,6)
    
    if ~exist('threshold', 'var')
        threshold = floor(max(max(image))/5);
    end
    if ~exist('particleLength', 'var')
        particleLength = 5.0;
    end
    if ~exist('filter', 'var')
        filter = 'bandpass';
    end
    if ~exist('noiseLength', 'var')
        noiseLength = 1.0;
    end
    if ~exist('visualize', 'var')
        visualize = ' ';
    end
    if(strcmp(filter,'bandpass'))
        %Gaussian filter
        G = gauss(10*noiseLength+1,10*noiseLength+1,noiseLength);

        %Box blur filter
        B = ones(round(particleLength), round(particleLength));
        B = B./sum(B(:));

        gConv = conv2(image,G,'same');
        bConv = conv2(image,B,'same');

        filtered = gConv - bConv;

        filtered(1:(round(particleLength)),:) = 0;
        filtered((end - particleLength + 1):end,:) = 0;
        filtered(:,1:(round(particleLength))) = 0;
        filtered(:,(end - particleLength + 1):end) = 0;
        filtered(filtered < 0) = 0;

    elseif(strcmp(filter, 'lowpass'))
        %Gaussian filter
        G = gauss(10*noiseLength+1,10*noiseLength+1,noiseLength);
        filtered = conv2(image,G,'same');
        filtered(1:(round(particleLength)),:) = 0;
        filtered((end - particleLength + 1):end,:) = 0;
        filtered(:,1:(round(particleLength))) = 0;
        filtered(:,(end - particleLength + 1):end) = 0;
        filtered(filtered < 0) = 0;
    else
        filtered = image;
    end

    %Finding all the pixels above threshold
    [row, col]=find(filtered > threshold);
    [nr,nc]=size(filtered);
    n=length(row);
    if n==0
        coordinates=[];
        disp('nothing above threshold');
        return;
    end

    rc=[row, col];
    brightPoints=zeros(size(row,1), 2);
    j = 0;
    for i=1:n
        r=rc(i,1);c=rc(i,2);
        %Checking each pixel above threshold to see if it's brighter than it's
        %neighbors.
        if r>1 && r<nr && c>1 && c<nc
            if (filtered(r,c)>=filtered(r-1,c-1) && filtered(r,c)>=filtered(r,c-1) ...
                    && filtered(r,c)>=filtered(r+1,c-1) && ...
                    filtered(r,c)>=filtered(r-1,c)  && filtered(r,c)>=filtered(r+1,c) ...
                    && filtered(r,c)>=filtered(r-1,c+1) && ...
                    filtered(r,c)>=filtered(r,c+1) && filtered(r,c)>=filtered(r+1,c+1))
                j = j + 1;
                brightPoints(j, 1) = c;
                brightPoints(j, 2) = r;
            end
        end
    end

    %Removing excess zeros.
    brightPoints(all(~brightPoints,2), : ) = [];

    [n,~]=size(brightPoints);
    %Getting rid of points that overlap.
    overlappingPoints = zeros(n,2);
    overlappingPointsIndexes = zeros(n,1);
    pointsToRemove = zeros(n,1);
    l = 2;
    %Looping over every bright point and comparing it to all the others.
    for i=1:n-1
        %This makes sure that every bright point is compared to another
        %bright point only once and that it doesn't get compared to itself.
        for j=i+1:n
            %Checking if two bright points overlap.
            if(((brightPoints(i,1)-brightPoints(j,1))^2 + ...
                    (brightPoints(i,2)-brightPoints(j,2))^2) < (2*particleLength)^2)
                %Saving the overlapping points and accompanying indexes.
                overlappingPoints(1,1) = brightPoints(i,1);
                overlappingPoints(1,2) = brightPoints(i,2);
                overlappingPoints(l,1) = brightPoints(j,1);
                overlappingPoints(l,2) = brightPoints(j,2);
                overlappingPointsIndexes(1) = i;
                overlappingPointsIndexes(l) = j;
                l = l + 1;
            end
        end
        %Removing excess zeros.
        overlappingPoints(all(~overlappingPoints,2), : ) = [];
        overlappingPointsIndexes(all(~overlappingPointsIndexes,2), : ) = [];
        l = 2;
        if(size(overlappingPoints,1) ~= 0)
            if(size(overlappingPoints,1) > 1)
                %Comparing overlapping points and keeping only the brightest
                %point.
                brightestPointIndex = 1;
                for k=1:size(overlappingPoints,1)-1
                    if(image(overlappingPoints(brightestPointIndex,1), ...
                            overlappingPoints(brightestPointIndex,2)) > ...
                            image(overlappingPoints(k+1,1),overlappingPoints(k+1,2)))
                    else
                        brightestPointIndex = k + 1;
                    end
                end
                %Making sure everything apart from the brightest point gets
                %removed.
                m = 1;
                for k=1:size(overlappingPoints,1)
                    if(k ~= brightestPointIndex)
                        pointsToRemove(m) = k;
                        m = m + 1;
                    end
                end
                %Removing excess zeros.
                pointsToRemove(all(~pointsToRemove,2), : ) = [];

                %Setting the less bright points to zero.
                for k=1:size(pointsToRemove,1)
                    brightPoints(overlappingPointsIndexes(pointsToRemove(k)),1) = 0;
                    brightPoints(overlappingPointsIndexes(pointsToRemove(k)),2) = 0;
                end
            end
        end
        %Resetting vectors.
        overlappingPoints = zeros(n,2);
        overlappingPointsIndexes = zeros(n,1);
        pointsToRemove = zeros(n,1);
    end
    %Removing the zero rows.
    brightPoints( ~any(brightPoints,2), : ) = [];

    %Plotting the coordinates onto the original given image.
    if(strcmp(visualize,'plot'))
        figure
        imagesc(image)
        hold on
        for i=1:1:size(brightPoints,1)
            viscircles([brightPoints(i,1) brightPoints(i,2)], 1,'Color', 'w', 'LineWidth', 2);
        end
        hold off
    end

    coordinates = brightPoints;
end