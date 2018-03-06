function [varargout]=PTV(images, masks, threshold, particleLength, filter, noiseLength, overlap, subWindow, searchArea, trackingSearchRadius, replaceOutliers, pivOnce, multipass)
    
    %Minimum 1 inputs and maximum 12 inputs
    narginchk(1,12)
    if ~exist('masks', 'var')
        masks = [];
    else
        %extending the masks to the amount of images.
        if(size(masks,2) < size(images,2))
            masks(end:numel(images))=masks(1,end);
        end
    end
    if ~exist('threshold', 'var')
        overlap = 0.5;
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
    if ~exist('overlap', 'var')
        overlap = 0.5;
    end
    if ~exist('subWindow', 'var')
        subWindow = floor(particleLength*12);
    end
    if ~exist('searchArea', 'var')
        searchArea = floor(subWindow/3);
    end
    if ~exist('trackingSearchRadius', 'var')
        trackingSearchRadius = particleLength;
    end
    if ~exist('replaceOutliers', 'var')
        replaceOutliers = 2;
    end
    if ~exist('pivOnce', 'var')
        pivOnce = 0;
    end
    if ~exist('multipass', 'var')
        multipass = [];
    end
    
    
    uniqueID = 1;
    firstFoundTrackingPoints = [];
    duplicateTrackingPoints = [];
    lastImageFFTP = [];
    lastImageDTP = [];
    firstTwoImages = 1;
    oneTime = 0;
    
    %For ensemble piv once.
    if(pivOnce == 1)
        %Setting up the loading dialog window
        waitBar = waitbar(0,'Running ensemble PIV...','Name','Particle Tracking Velocimetry.');
        
        im1 = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
        im2 = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
        n = 1;
        for imageNumber=1:2:size(images,2)-1
            im1(:,:,n) = im2uint8(imread(images{imageNumber}));
            im2(:,:,n) = im2uint8(imread(images{imageNumber + 1}));
            n = n + 1;
        end
        opt = setpivopt('range',[-searchArea searchArea -searchArea searchArea],'subwindow',subWindow,subWindow,overlap,'ensemble',@nanmedian);
        if isempty(masks)
            piv = normalpass(multipass,im1,masks,im2,masks,opt);
        else
            mask1 = imread(masks{1});
            mask2 = imread(masks{1});
            piv = normalpass(multipass,im1,mask1,im2,mask1,opt);
        end
        [U,V] = replaceoutliers(piv);
        oneTime = 1;
    else
        %Setting up the loading dialog window
        waitBar = waitbar(0,'Please wait... PTV 0.0% finished.','Name','Particle Tracking Velocimetry.');
    end
    
    
    for imageNumber=1:size(images,2)-1
        loadingPercentText = ['Please wait... PTV ', sprintf('%.1f',imageNumber*100.0 / (size(images,2))), '% finished.'];
        waitbar(imageNumber / (size(images,2)), waitBar, loadingPercentText)
       
        im1 = im2uint8(imread(images{imageNumber}));
        im2 = im2uint8(imread(images{imageNumber + 1}));
        if(pivOnce == 0)
            if isempty(masks)
                [U, V] = PIV(im1, im2, masks, masks, overlap, subWindow, searchArea, replaceOutliers, multipass);
            else
                mask1 = imread(masks{imageNumber});
                mask2 = imread(masks{imageNumber + 1});
                [U, V] = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, replaceOutliers, multipass);
            end
        end
        if(firstTwoImages == 1)
            coordinates1 = findParticles(im1, threshold, particleLength, filter, noiseLength);
            coordinates2 = findParticles(im2, threshold, particleLength, filter, noiseLength);
        else
            coordinates1 = coordinates2;
            coordinates2 = findParticles(im2, threshold, particleLength, filter, noiseLength);
        end

        %Removing particles outside of mask
        if(firstTwoImages == 1)
            if ~isempty(mask1) && ~isempty(mask2)
                coordinatesMasked1 = coordinates1;
                indexToDelete = [];
                for i=1:size(coordinatesMasked1,1)
                    if (mask1(coordinates1(i,2),coordinates1(i,1)) == 0)
                        indexToDelete = [indexToDelete, i];
                    end
                end
                if ~isempty(indexToDelete)
                    coordinatesMasked1(indexToDelete,:) = [];
                end

                coordinatesMasked2 = coordinates2;
                indexToDelete = [];
                for i=1:size(coordinatesMasked2,1)
                    if (mask2(coordinates2(i,2),coordinates2(i,1)) == 0)
                        indexToDelete = [indexToDelete, i];
                    end
                end
                if ~isempty(indexToDelete)
                    coordinatesMasked2(indexToDelete,:) = [];
                end
            end
        else
            coordinatesMasked1 = coordinatesMasked2;
            if ~isempty(mask1) && ~isempty(mask2)
                coordinatesMasked2 = coordinates2;
                indexToDelete = [];
                for i=1:size(coordinatesMasked2,1)
                    if (mask2(coordinates2(i,2),coordinates2(i,1)) == 0)
                        indexToDelete = [indexToDelete, i];
                    end
                end
                if ~isempty(indexToDelete)
                    coordinatesMasked2(indexToDelete,:) = [];
                end
            end
        end
        
        if(pivOnce == 0)
            %Interpolating the vertical and horizontal velocities up to the size of
            %the images
            x = linspace(0, size(U,1), size(U,1));
            y = linspace(0, size(U,2), size(U,2));
            xq = linspace(0, size(U,1), size(im1,1));
            yq = linspace(0, size(U,2), size(im1,2));
            [X, Y] = meshgrid(x,y);
            [Xq, Yq] = meshgrid(xq,yq);
            Uint = interp2(X,Y,U,Xq,Yq, 'cubic');

            x = linspace(0, size(V,1), size(V,1));
            y = linspace(0, size(V,2), size(V,2));
            xq = linspace(0, size(V,1), size(im1,1));
            yq = linspace(0, size(V,2), size(im1,2));
            [X, Y] = meshgrid(x,y);
            [Xq, Yq] = meshgrid(xq,yq);
            Vint = interp2(X,Y,V,Xq,Yq, 'cubic');
        elseif(oneTime == 1)
            %Interpolating the vertical and horizontal velocities up to the size of
            %the images
            x = linspace(0, size(U,1), size(U,1));
            y = linspace(0, size(U,2), size(U,2));
            xq = linspace(0, size(U,1), size(im1,1));
            yq = linspace(0, size(U,2), size(im1,2));
            [X, Y] = meshgrid(x,y);
            [Xq, Yq] = meshgrid(xq,yq);
            Uint = interp2(X,Y,U,Xq,Yq, 'cubic');

            x = linspace(0, size(V,1), size(V,1));
            y = linspace(0, size(V,2), size(V,2));
            xq = linspace(0, size(V,1), size(im1,1));
            yq = linspace(0, size(V,2), size(im1,2));
            [X, Y] = meshgrid(x,y);
            [Xq, Yq] = meshgrid(xq,yq);
            Vint = interp2(X,Y,V,Xq,Yq, 'cubic');
            oneTime = 0;
        end
        %Finding new particle positions from vertical and horizontal
        %velocity fields.
        xInt = zeros(size(coordinatesMasked1,1),1);
        yInt = zeros(size(coordinatesMasked1,1),1);
        for i=1:size(coordinatesMasked1,1)
            xInt(i) = Uint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)) + coordinatesMasked1(i,1);
            yInt(i) = Vint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)) + coordinatesMasked1(i,2);
        end
        
        %Connecting particles from different images together.
        firstTime = 1;
        if(firstTwoImages == 1)
            %Adding the coordinates of both images for the two 
            %first images.
            for i=1:size(coordinatesMasked1,1)
                for j=1:size(coordinatesMasked2,1)
                    %Checking if a interpolated particle is within its own
                    %radius in the second image. If yes then we will connect
                    %the particle from the first image to the particle in the
                    %second image.
                    %IMPROVEMENTS TO BE MADE: Could make the radius be moved in
                    %the direction of the flow so that we get a higher
                    %probability of finding the correct tracking.
                    if((xInt(i) - coordinatesMasked2(j,1))^2 + (yInt(i) - coordinatesMasked2(j,2))^2 < trackingSearchRadius^2)
                        if(firstTime == 1)
                            firstTime = 0;
                            firstFoundTrackingPoints = [firstFoundTrackingPoints; xInt(i) - Uint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)), yInt(i) - Vint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)), imageNumber, uniqueID];
                            firstFoundTrackingPoints = [firstFoundTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                            %Saving the coordinates of the most recent
                            %image in a seperate list for use in the next
                            %for loop increment.
                            lastImageFFTP = [lastImageFFTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                            uniqueID = uniqueID + 1;
                        else
                            %More than one particle may be found within that
                            %radius. We will only keep the first one. Hopefully
                            %this list is empty.
                            duplicateTrackingPoints = [duplicateTrackingPoints; xInt(i) - Uint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)), yInt(i) - Vint(coordinatesMasked1(i,2), coordinatesMasked1(i,1)), imageNumber, uniqueID];
                            duplicateTrackingPoints = [duplicateTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                            %Saving the coordinates of the most recent
                            %image in a seperate list for use in the next
                            %for loop increment.
                            lastImageDTP = [lastImageDTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                        end
                    end
                end
                firstTime = 1;
            end
            firstTwoImages = 0;
        else
            staticLIFFTP = lastImageFFTP;
            staticLIDTP = lastImageDTP;
            %Resetting them so that they can be filled up with the information
            %of the next latest image.
            lastImageFFTP = [];
            lastImageDTP = [];
            %This is for keeping track of if we have found any particles in
            %our two above lists.
            particleInArray = 0;
            for i=1:size(coordinatesMasked1,1)
                for j=1:size(coordinatesMasked2,1)
                    %Checking if a interpolated particle is within its own
                    %radius in the second image. If yes then we will connect
                    %the particle from the first image to the particle in the
                    %second image.
                    %IMPROVEMENTS TO BE MADE: Could make the radius be moved in
                    %the direction of the flow so that we get a higher
                    %probability of finding the correct tracking.
                    if((xInt(i) - coordinatesMasked2(j,1))^2 + (yInt(i) - coordinatesMasked2(j,2))^2 < trackingSearchRadius^2)
                        for k=1:size(staticLIFFTP,1)
                            if((staticLIFFTP(k,1) - (xInt(i) - Uint(coordinatesMasked1(i,2), coordinatesMasked1(i,1))))^2 + (staticLIFFTP(k,2) - (yInt(i) - Vint(coordinatesMasked1(i,2), coordinatesMasked1(i,1))))^2 < (trackingSearchRadius/4)^2)
                                if(firstTime == 1)
                                    particleInArray = 1;
                                    firstTime = 0;
                                    firstFoundTrackingPoints = [firstFoundTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                    lastImageFFTP = [lastImageFFTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                else
                                    particleInArray = 1;
                                    %More than one particle may be found within that
                                    %radius. We will only keep the first one. Hopefully
                                    %this list is empty.
                                    duplicateTrackingPoints = [duplicateTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                    lastImageDTP = [lastImageDTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                end
                            end
                        end
                        if(particleInArray == 0)
                            if(~isempty(staticLIDTP))
                                for l=1:size(staticLIDTP,1)
                                    if((staticLIDTP(l,1) - (xInt(i) - Uint(coordinatesMasked1(i,2), coordinatesMasked1(i,1))))^2 + (staticLIDTP(l,2) - (yInt(i) - Vint(coordinatesMasked1(i,2), coordinatesMasked1(i,1))))^2 < (trackingSearchRadius/4)^2)
                                        if(firstTime == 1)
                                            particleInArray = 1;
                                            firstTime = 0;
                                            firstFoundTrackingPoints = [firstFoundTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                            lastImageFFTP = [lastImageFFTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                        else
                                            particleInArray = 1;
                                            %More than one particle may be found within that
                                            %radius. We will only keep the first one. Hopefully
                                            %this list is empty.
                                            duplicateTrackingPoints = [duplicateTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                            lastImageDTP = [lastImageDTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, staticLIFFTP(k,4)];
                                        end
                                    end
                                end
                            end
                        end
                        if(particleInArray == 0)
                            if(firstTime == 1)
                                particleInArray = 1;
                                firstTime = 0;
                                firstFoundTrackingPoints = [firstFoundTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                                %Saving the coordinates of the most recent
                                %image in a seperate list for use in the next
                                %for loop increment.
                                lastImageFFTP = [lastImageFFTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                                uniqueID = uniqueID + 1;
                            else
                                particleInArray = 1;
                                %More than one particle may be found within that
                                %radius. We will only keep the first one. Hopefully
                                %this list is empty.
                                duplicateTrackingPoints = [duplicateTrackingPoints; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                                %Saving the coordinates of the most recent
                                %image in a seperate list for use in the next
                                %for loop increment.
                                lastImageDTP = [lastImageDTP; coordinatesMasked2(j,1), coordinatesMasked2(j,2), imageNumber+1, uniqueID];
                            end
                        end
                    end
                    particleInArray = 0;
                end
                firstTime = 1;
            end
        end
    end
    
    close(waitBar)
    msgbox('Operation Completed');
    
    %return
    if nargout == 1
        varargout{1} = firstFoundTrackingPoints;
    elseif nargout == 2
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = duplicateTrackingPoints;
    end
end