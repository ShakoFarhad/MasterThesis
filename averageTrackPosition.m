function averageTP = averageTrackPosition(tracks)
    IDsortedTracks = sortrows(tracks,4);
    findDifferenceIndex = find(diff(IDsortedTracks(:,4)) ~= 0);
    trajectoryLengthsInArray = diff(findDifferenceIndex);
    
    averageTP = zeros(size(trajectoryLengthsInArray,1)+2,2);
    
    x = mean(IDsortedTracks(1:findDifferenceIndex(1),1));
    y = mean(IDsortedTracks(1:findDifferenceIndex(1),2));
    averageTP(1,:) = [x, y];

    for i=1:size(trajectoryLengthsInArray,1)
        x = mean(IDsortedTracks(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),1));
        y = mean(IDsortedTracks(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),2));
        averageTP(i+1,:) = [x, y];
    end
    x = mean(IDsortedTracks(findDifferenceIndex(end)+1:size(IDsortedTracks,1),1));
    y = mean(IDsortedTracks(findDifferenceIndex(end)+1:size(IDsortedTracks,1),2));
    averageTP(end,:) = [x, y];
end