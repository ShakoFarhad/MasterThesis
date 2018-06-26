function averageTV = averageTrackVelocity(tracks)
    IDsortedTracks = sortrows(tracks,4);
    findDifferenceIndex = find(diff(IDsortedTracks(:,4)) ~= 0);
    trajectoryLengthsInArray = diff(findDifferenceIndex);
    
    averageTV = zeros(size(trajectoryLengthsInArray,1)+2,2);
    
    x = mean(IDsortedTracks(1:findDifferenceIndex(1),6));
    y = mean(IDsortedTracks(1:findDifferenceIndex(1),7));
    
    averageTV(1,:) = [x, y];
    for i=1:size(trajectoryLengthsInArray,1)
        x = mean(IDsortedTracks(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),6));
        y = mean(IDsortedTracks(findDifferenceIndex(i)+1:findDifferenceIndex(i+1),7));
        averageTV(i+1,:) = [x, y];
    end
    x = mean(IDsortedTracks(findDifferenceIndex(end)+1:size(IDsortedTracks,1),6));
    y = mean(IDsortedTracks(findDifferenceIndex(end)+1:size(IDsortedTracks,1),7));
    averageTV(end,:) = [x, y];
end