function [varargout]=PTV(images, masks, threshold, particleLength, filter, noiseLength, overlap, subWindow, searchArea, distortedPass, trackingSearchRadius, particleSkippedFrames, replaceOutliers, ensembleContribution, simplifyContribution, previousVelocityContribution, blur, sigma, subpixel)
    
    %Minimum 1 inputs and maximum 19 inputs
    narginchk(1,19)
    if ~exist('masks', 'var')
        masks = [];
    elseif ~isempty(masks)
        %extending the masks to the amount of images.
        if(size(masks,2) < size(images,2))
            masks(end:numel(images))=masks(1,end);
        end
    end
    if ~exist('threshold', 'var')
        threshold = mean(mean(imread(images{1})));
    end
    if ~exist('particleLength', 'var')
        particleLength = 5.0;
    end
    if ~exist('filter', 'var')
        filter = 'Aggressive';
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
        searchArea = floor(min(subWindow)/3);
    end
    if ~exist('distortedPass', 'var')
        distortedPass = 0;
    end
    if ~exist('trackingSearchRadius', 'var')
        trackingSearchRadius = 2;
    end
    if ~exist('particleSkippedFrames', 'var')
        particleSkippedFrames = 0;
    end
    if ~exist('replaceOutliers', 'var')
        replaceOutliers = 2;
    end
    if ~exist('ensembleContribution', 'var')
        ensembleContribution = 0;
    else
        if(ensembleContribution > 1)
            ensembleContribution = 1;
        elseif(ensembleContribution < 0)
            ensembleContribution = 0;
        end
    end
    if ~exist('simplifyContribution', 'var')
        simplifyContribution = 0;
    else
        if(simplifyContribution > 1)
            simplifyContribution = 1;
        elseif(simplifyContribution < 0)
            simplifyContribution = 0;
        end
    end
    if ~exist('previousVelocityContribution', 'var')
        previousVelocityContribution = 0;
    else
        if(previousVelocityContribution > 1)
            previousVelocityContribution = 1;
        elseif(previousVelocityContribution < 0)
            previousVelocityContribution = 0;
        end
    end
    if ~exist('blur', 'var')
        blur = 'none';
    end
    if ~exist('sigma', 'var')
        sigma = 0.5;
    end
    if ~exist('subpixel', 'var')
        subpixel = @subpixelnone;
    end
    
    multipass = [];
    
    uniqueID = 1;
    %Finding out how big the arrays can be. And allocating 1% of
    %total/available memory to these arrays.
    if(ispc)
        user = memory;
        doubleArraySize = round((user.MaxPossibleArrayBytes/8)/200);
        firstFoundTrackingPoints = zeros(doubleArraySize,7);
        lastImageFFTP = zeros(doubleArraySize,7);
    else
        firstFoundTrackingPoints = zeros(2^21,7);
        lastImageFFTP = zeros(2^21,7);
    end
    FFTPCounter = 1;
    LIFFTPCounter = 1;
    firstTwoImages = 1;
    
    %These values are used to keep track of skipped particles, or
    %particles that dissappear for N frames, (particlesSkippedFrames).
    keepParticleTrack = 0;
    guessedParticleTrack = 1;
    
    %Squaring the search radius so that we do not have to do it for every
    %iteration to save computational power.
    trackingSearchRadiusSquared = trackingSearchRadius^2;
    
    %For ensemble piv once.
    if(ensembleContribution > 0 && ensembleContribution <= 1)
        %Setting up the loading dialog window
        waitBar = waitbar(0,'Running ensemble PIV...','Name','Particle Tracking Velocimetry.');
        
        im1PIV = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
        im2PIV = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
        n = 1;
        for imageNumber=1:2:size(images,2)-1
            tempImage1 = imread(images{imageNumber});
            tempImage2 = imread(images{imageNumber + 1});
            if(size(tempImage1,3) == 3)
                im1 = im2double(rgb2gray(tempImage1));
            else
                im1 = im2double(tempImage1(:,:,1));
            end
            if(size(tempImage2,3) == 3)
                im2 = im2double(rgb2gray(tempImage2));
            else
                im2 = im2double(tempImage2(:,:,1));
            end
            coordinates1 = findParticles(im1, threshold, particleLength, filter, noiseLength);
            coordinates2 = findParticles(im2, threshold, particleLength, filter, noiseLength);
            %Removing particles outside of mask
            if ~isempty(masks)
                mask1 = imread(masks{1});
                mask2 = imread(masks{1});
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
            if(simplifyContribution == 1)
                drawIm1 = im2double(drawCircles(im1, coordinates1, particleLength, blur, sigma));
                drawIm2 = im2double(drawCircles(im2, coordinates2, particleLength, blur, sigma));
                im1PIV(:,:,n) = drawIm1;
                im2PIV(:,:,n) = drawIm2;
            elseif(simplifyContribution == 0)
                im1PIV(:,:,n) = im1;
                im2PIV(:,:,n) = im2;
            else
                drawIm1 = im2double(drawCircles(im1, coordinates1, particleLength, blur, sigma));
                drawIm2 = im2double(drawCircles(im2, coordinates2, particleLength, blur, sigma));
                im1PIV(:,:,n) = im1*(1-simplifyContribution) + drawIm1/max(max(drawIm1))*max(max(im1))*simplifyContribution;
                im2PIV(:,:,n) = im2*(1-simplifyContribution) + drawIm2/max(max(drawIm2))*max(max(im2))*simplifyContribution;
            end
            n = n + 1;
        end
        if(size(searchArea,2) == 1)
            x1 = searchArea(1);
            y1 = searchArea(1);
        else
            x1 = searchArea(1);
            y1 = searchArea(2);
        end
        if(size(subWindow,2) == 1)
            x = subWindow(1);
            y = subWindow(1);
        else
            x = subWindow(1);
            y = subWindow(2);
        end
        opt = setpivopt('range',[-x1 x1 -y1 y1],'subwindow',x,y,overlap,'ensemble',@nanmedian, 'subpixel', subpixel);
        if isempty(masks)
            piv = normalpass(multipass,im1PIV,masks,im2PIV,masks,opt);
            if(distortedPass>0)
                %Shifted pass using smaller subwindows
                multipass = piv;
                piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 0, multipass, subpixel, particleLength);
                % Distorted passes
                for distortedPasses=1:distortedPass-1
                    multipass = piv;
                    piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 1, multipass, subpixel, particleLength);
                end
            end
        else
            mask1 = imread(masks{1});
            piv = normalpass(multipass,im1PIV,mask1,im2PIV,mask1,opt);
            if(distortedPass>0)
                %Shifted pass using smaller subwindows
                multipass = piv;
                piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 0, multipass, subpixel, particleLength);
                % Distorted passes
                for distortedPasses=1:distortedPass-1
                    multipass = piv;
                    piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 1, multipass, subpixel, particleLength);
                end
            end
        end
        if(replaceOutliers == 1)
            [U,V] = replaceoutliers(piv,mask1&mask1);
            piv.U = U;
            piv.V = V;
        elseif(replaceOutliers == 2)
            [U,V] = replaceoutliers(piv);
            piv.U = U;
            piv.V = V;
        else
            U = piv.U;
            V = piv.V;
        end
        x = linspace(0, size(U,1), size(U,1));
        y = linspace(0, size(U,2), size(U,2));
        xq = linspace(0, size(U,1), size(im1,1));
        yq = linspace(0, size(U,2), size(im1,2));
        [X, Y] = meshgrid(y,x);
        [Xq, Yq] = meshgrid(yq,xq);
        ensembleUint = interp2(X,Y,U,Xq,Yq, 'cubic');

        x = linspace(0, size(V,1), size(V,1));
        y = linspace(0, size(V,2), size(V,2));
        xq = linspace(0, size(V,1), size(im1,1));
        yq = linspace(0, size(V,2), size(im1,2));
        [X, Y] = meshgrid(y,x);
        [Xq, Yq] = meshgrid(yq,xq);
        ensembleVint = interp2(X,Y,V,Xq,Yq, 'cubic');
        if(ensembleContribution == 1)
            %Interpolating the velocity fields only once
            Uint = ensembleUint;
            Vint = ensembleVint;
        end
    else
        ensembleUint = 0; ensembleVint = 0;
        %Setting up the loading dialog window
        waitBar = waitbar(0,'Please wait... PTV 0.0% finished.','Name','Particle Tracking Velocimetry.');
    end
    
    
    for imageNumber=1:size(images,2)-1
        try
            loadingPercentText = ['Please wait... PTV ', sprintf('%.1f',imageNumber*100.0 / (size(images,2))), '% finished.'];
            waitbar(imageNumber / (size(images,2)), waitBar, loadingPercentText)
        catch
            disp("Loading bar cancelled. Exiting PTV.");
            break
        end
       
        tempImage1 = imread(images{imageNumber});
        tempImage2 = imread(images{imageNumber + 1});
        if(size(tempImage1,3) == 3)
            im1 = im2double(rgb2gray(tempImage1));
        else
            im1 = im2double(tempImage1(:,:,1));
        end
        if(size(tempImage1,3) == 3)
            im2 = im2double(rgb2gray(tempImage2));
        else
            im2 = im2double(tempImage2(:,:,1));
        end
        
        if(firstTwoImages == 1)
            coordinates1 = findParticles(im1, threshold, particleLength, filter, noiseLength);
            coordinates2 = findParticles(im2, threshold, particleLength, filter, noiseLength);
            coordinates1 = [coordinates1 zeros(size(coordinates1,1),1)];
            coordinates2 = [coordinates2 zeros(size(coordinates2,1),1)];
        else
            coordinates1 = coordinates2;
            coordinates2 = findParticles(im2, threshold, particleLength, filter, noiseLength);
            coordinates2 = [coordinates2 zeros(size(coordinates2,1),1)];
        end
        
        if ~isempty(masks)
            mask1 = imread(masks{imageNumber});
            mask2 = imread(masks{imageNumber + 1});
        else
            mask1 = masks;
            mask2 = masks;
        end

        %Removing particles outside of mask
        if(firstTwoImages == 1)
            if(~isempty(mask1) && ~isempty(mask2))
                coordinatesMasked1 = coordinates1;
                indexToDelete = [];
                for i=1:size(coordinatesMasked1,1)
                    if(mask1(coordinates1(i,2),coordinates1(i,1)) == 0)
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
            if ~isempty(mask1) && ~isempty(mask2)
                coordinatesMasked1 = coordinatesMasked2;
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
        if isempty(masks)
            coord1 = coordinates1;
            coord2 = coordinates2;
        else
            coord1 = coordinatesMasked1;
            coord2 = coordinatesMasked2;
        end

        %Simplifying the images for the PIV by drawing circles where particles have been found.    
        if(simplifyContribution == 1)
            im1 = im2double(drawCircles(im1, coord1(:,1:2), particleLength, blur, sigma));
            im2 = im2double(drawCircles(im2, coord2(:,1:2), particleLength, blur, sigma));
            if(distortedPass>0)
                multipass = [];
                piv = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, 0, 0, multipass,subpixel, particleLength);
                %Shifted pass using smaller subwindows
                multipass = piv;
                piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 0, multipass,subpixel, particleLength);
                % Distorted passes
                for distortedPasses=1:distortedPass-1
                    multipass = piv;
                    piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 1, multipass,subpixel, particleLength);
                end
                multipass = piv;
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), replaceOutliers, 1, multipass,subpixel, particleLength);
            else
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, replaceOutliers, 1, multipass,subpixel, particleLength);
            end
        elseif(simplifyContribution == 0)
            if(distortedPass>0)
                multipass = [];
                piv = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, 0, 0, multipass,subpixel, particleLength);
                %Shifted pass using smaller subwindows
                multipass = piv;
                piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 0, multipass,subpixel, particleLength);
                % Distorted passes
                for distortedPasses=1:distortedPass-1
                    multipass = piv;
                    piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 1, multipass,subpixel, particleLength);
                end
                multipass = piv;
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), replaceOutliers, 1, multipass,subpixel, particleLength);
            else
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, replaceOutliers, distortedPass, multipass,subpixel, particleLength);
            end
        else
            drawIm1 = im2double(drawCircles(im1, coord1(:,1:2), particleLength, blur, sigma));
            drawIm2 = im2double(drawCircles(im2, coord2(:,1:2), particleLength, blur, sigma));
            im1 = im1*(1-simplifyContribution) + drawIm1/max(max(drawIm1))*max(max(im1))*simplifyContribution;
            im2 = im2*(1-simplifyContribution) + drawIm2/max(max(drawIm2))*max(max(im2))*simplifyContribution;
            if(distortedPass>0)
                multipass = [];
                piv = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, 0, 0, multipass,subpixel, particleLength);
                %Shifted pass using smaller subwindows
                multipass = piv;
                piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 0, multipass,subpixel, particleLength);
                % Distorted passes
                for distortedPasses=1:distortedPass-1
                    multipass = piv;
                    piv = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), 0, 1, multipass,subpixel, particleLength);
                end
                multipass = piv;
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, round(subWindow/2), round(searchArea/2), replaceOutliers, 1, multipass,subpixel, particleLength);
            else
                [piv, U, V] = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, replaceOutliers, distortedPass, multipass,subpixel, particleLength);
            end
        end
        
        if(ensembleContribution >= 0 && ensembleContribution < 1)
            %Interpolating the vertical and horizontal velocities up to the size of
            %the images
            x = linspace(0, size(U,1), size(U,1));
            y = linspace(0, size(U,2), size(U,2));
            xq = linspace(0, size(U,1), size(im1,1));
            yq = linspace(0, size(U,2), size(im1,2));
            [X, Y] = meshgrid(y,x);
            [Xq, Yq] = meshgrid(yq,xq);
            Uint = interp2(X,Y,U,Xq,Yq, 'cubic');

            x = linspace(0, size(V,1), size(V,1));
            y = linspace(0, size(V,2), size(V,2));
            xq = linspace(0, size(V,1), size(im1,1));
            yq = linspace(0, size(V,2), size(im1,2));
            [X, Y] = meshgrid(y,x);
            [Xq, Yq] = meshgrid(yq,xq);
            Vint = interp2(X,Y,V,Xq,Yq, 'cubic');
            if(ensembleContribution ~= 0)
                Uint = Uint*(1-ensembleContribution) + ensembleUint*ensembleContribution;
                Vint = Vint*(1-ensembleContribution) + ensembleVint*ensembleContribution;
            end
        end
        
        %Connecting particles from different images together.
        if(firstTwoImages == 1)
            %Finding new particle positions from vertical and horizontal
            %velocity fields.
            xInt = zeros(size(coord1,1),1);
            yInt = zeros(size(coord1,1),1);
            for i=1:size(coord1,1)
                %We do not want to save points that do not move at all.
                if(~(abs(Uint(coord1(i,2), coord1(i,1))) < 0.5 && abs(Vint(coord1(i,2), coord1(i,1))) < 0.5))
                    xInt(i) = Uint(coord1(i,2), coord1(i,1)) + coord1(i,1);
                    yInt(i) = Vint(coord1(i,2), coord1(i,1)) + coord1(i,2);
                end
            end
            
            %Adding the coordinates of both images for the two 
            %first images.
            for i=1:size(coord1,1)
                for j=1:size(coord2,1)
                    %Checking if a interpolated particle is within the
                    %search radius in the second image. If yes then we will
                    %connect the particle from the first image to the 
                    %particle in the second image.
                    if(xInt(i) ~= 0)
                        if((xInt(i) - coord2(j,1))^2 + (yInt(i) - coord2(j,2))^2 <= trackingSearchRadiusSquared)
                            coord1(i,3) = 1; coord2(j,3) = 1;
                            firstFoundTrackingPoints(FFTPCounter,:) = [coord1(i,1), coord1(i,2), imageNumber, uniqueID, keepParticleTrack, Uint(coord1(i,2), coord1(i,1)), Vint(coord1(i,2), coord1(i,1))];
                            FFTPCounter = FFTPCounter + 1;
                            firstFoundTrackingPoints(FFTPCounter,:) = [coord2(j,1), coord2(j,2), imageNumber+1, uniqueID, keepParticleTrack, Uint(coord1(i,2), coord1(i,1)), Vint(coord1(i,2), coord1(i,1))];
                            FFTPCounter = FFTPCounter + 1;
                            %Saving the coordinates of the most recent
                            %image in a seperate list for use in the next
                            %image for loop increment.
                            lastImageFFTP(LIFFTPCounter,:) = [coord2(j,1), coord2(j,2), imageNumber+1, uniqueID, keepParticleTrack, Uint(coord1(i,2), coord1(i,1)), Vint(coord1(i,2), coord1(i,1))];
                            LIFFTPCounter = LIFFTPCounter + 1;
                            uniqueID = uniqueID + 1;
                            break
                        end
                    end
                end
            end
            %Extracting all the coordinates that were not connected in a
            %track.
            coord1 = coord1(coord1(:,3) == 0,:);
            coord2 = coord2(coord2(:,3) == 0,:);
            %We are guessing where the hidden particle will be in the next 
            %image.
            if(particleSkippedFrames > 0)
                for i=1:size(coord1,1)
                    if(xInt(i,1) ~= 0)
                        if(~isempty(mask1) && ~isempty(mask2))
                            %Making sure the interpolated particle is within the
                            %image and mask.
                            if(Uint(coord1(i,2), coord1(i,1)) + coord1(i,1) <= size(im1,2) && Uint(coord1(i,2), coord1(i,1)) + coord1(i,1) >= 1 && Vint(coord1(i,2), coord1(i,1)) + coord1(i,2) <= size(im1,1) && Vint(coord1(i,2), coord1(i,1)) + coord1(i,2) >= 1)
                                if(mask1(round(Vint(coord1(i,2), coord1(i,1)) + coord1(i,2)), round(Uint(coord1(i,2), coord1(i,1)) + coord1(i,1))) ~= 0)
                                    %Saving the coordinates of the most recent
                                    %image in a seperate list for use in the next
                                    %image for loop increment.
                                    lastImageFFTP(LIFFTPCounter, :) = [round(Uint(coord1(i,2), coord1(i,1)) + coord1(i,1)), round(Vint(coord1(i,2), coord1(i,1)) + coord1(i,2)), imageNumber+1, uniqueID, guessedParticleTrack, Uint(coord1(i,2), coord1(i,1)), Vint(coord1(i,2), coord1(i,1))];
                                    LIFFTPCounter = LIFFTPCounter + 1;
                                    uniqueID = uniqueID + 1;
                                end
                            end
                        else
                            if(Uint(coord1(i,2), coord1(i,1)) + coord1(i,1) <= size(im1,2) && Uint(coord1(i,2), coord1(i,1)) + coord1(i,1) >= 1 && Vint(coord1(i,2), coord1(i,1)) + coord1(i,2) <= size(im1,1) && Vint(coord1(i,2), coord1(i,1)) + coord1(i,2) >= 1)
                                %Saving the coordinates of the most recent
                                %image in a seperate list for use in the next
                                %image for loop increment.
                                lastImageFFTP(LIFFTPCounter, :) = [round(Uint(coord1(i,2), coord1(i,1)) + coord1(i,1)), round(Vint(coord1(i,2), coord1(i,1)) + coord1(i,2)), imageNumber+1, uniqueID, guessedParticleTrack, Uint(coord1(i,2), coord1(i,1)), Vint(coord1(i,2), coord1(i,1))];
                                LIFFTPCounter = LIFFTPCounter + 1;
                                uniqueID = uniqueID + 1;
                            end
                        end
                    end
                end
            end
            %We are adding all the new particles that could not be
            %connected to particles in the previous image.
            for j=1:size(coord2,1)
                %Saving the coordinates of the most recent
                %image in a seperate list for use in the next
                %image for loop increment.
                lastImageFFTP(LIFFTPCounter, :) = [coord2(j,1), coord2(j,2), imageNumber+1, uniqueID, keepParticleTrack, Uint(coord2(j,2), coord2(j,1)), Vint(coord2(j,2), coord2(j,1))];
                LIFFTPCounter = LIFFTPCounter + 1;
                uniqueID = uniqueID + 1;
            end
            firstTwoImages = 0;
        else
            staticLIFFTP = lastImageFFTP;
            %Resetting it so that it can be filled up with the information
            %of the next latest image.
            lastImageFFTP = zeros(size(staticLIFFTP,1)+size(coord2,1),7);
            LIFFTPCounter = 1;
            %Finding new particle positions from vertical and horizontal
            %velocity fields.
            lastImageInt = [staticLIFFTP zeros(size(staticLIFFTP,1),3)];
            for i=1:size(staticLIFFTP,1)
                prevCurrentVelocityU = Uint(staticLIFFTP(i,2), staticLIFFTP(i,1))*(1-previousVelocityContribution) + staticLIFFTP(i,6)*previousVelocityContribution;
                prevCurrentVelocityV = Vint(staticLIFFTP(i,2), staticLIFFTP(i,1))*(1-previousVelocityContribution) + staticLIFFTP(i,7)*previousVelocityContribution;
                lastImageInt(i,7) = prevCurrentVelocityU;
                lastImageInt(i,8) = prevCurrentVelocityV;
                lastImageInt(i,1) = prevCurrentVelocityU + staticLIFFTP(i,1);
                lastImageInt(i,2) = prevCurrentVelocityV + staticLIFFTP(i,2);
            end
            for i=1:size(lastImageInt,1)
                for j=1:size(coord2,1)
                    %Checking if a interpolated particle is within the
                    %search radius in the second image. If yes then we will
                    %connect the particle from the first image to the 
                    %particle in the second image.
                    if((lastImageInt(i,1) - coord2(j,1))^2 + (lastImageInt(i,2) - coord2(j,2))^2 <= trackingSearchRadiusSquared)
                        %Marking this particle as connected
                        lastImageInt(i,6) = 1; coord2(j,3) = 1;

                        firstFoundTrackingPoints(FFTPCounter, :) = [staticLIFFTP(i,1), staticLIFFTP(i,2), imageNumber, lastImageInt(i,4), keepParticleTrack, lastImageInt(i,7), lastImageInt(i,8)];
                        FFTPCounter = FFTPCounter + 1;
                        firstFoundTrackingPoints(FFTPCounter, :) = [coord2(j,1), coord2(j,2), imageNumber+1, lastImageInt(i,4), keepParticleTrack, Uint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Uint(coord2(j,2), coord2(j,1))*previousVelocityContribution, Vint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Vint(coord2(j,2), coord2(j,1))*previousVelocityContribution];
                        FFTPCounter = FFTPCounter + 1;

                        lastImageFFTP(LIFFTPCounter, :) = [coord2(j,1), coord2(j,2), imageNumber+1, lastImageInt(i,4), keepParticleTrack, Uint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Uint(coord2(j,2), coord2(j,1))*previousVelocityContribution, Vint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Vint(coord2(j,2), coord2(j,1))*previousVelocityContribution];
                        LIFFTPCounter = LIFFTPCounter + 1;
                        %We are only saving/connecting the first particle found within
                        %the radius.
                        break
                    end
                end
            end
            %Extracting all the coordinates that were not connected in a
            %track.
            lastImageInt = lastImageInt(lastImageInt(:,6) == 0,:);
            for i=1:size(lastImageInt,1)
                %We are guessing where the hidden particle will be in the next 
                %image.
                if(particleSkippedFrames > 0)
                    %Making sure the interpolated particle is within the
                    %image and mask.
                    if(~isempty(mask1) && ~isempty(mask2))
                        if(lastImageInt(i,1) <= size(im1,2) && lastImageInt(i,1) >= 1 && lastImageInt(i,2) <= size(im1,1) && lastImageInt(i,2) >= 1)
                            if(mask1(round(lastImageInt(i,2)), round(lastImageInt(i,1))) ~= 0)
                                firstFoundTrackingPoints(FFTPCounter, :) = [round(lastImageInt(i,1)), round(lastImageInt(i,2)), imageNumber, lastImageInt(i,4), guessedParticleTrack, lastImageInt(i,7), lastImageInt(i,8)];
                                FFTPCounter = FFTPCounter + 1;
                                lastImageFFTP(LIFFTPCounter, :) = [round(lastImageInt(i,1)), round(lastImageInt(i,2)), imageNumber+1, lastImageInt(i,4), guessedParticleTrack, lastImageInt(i,7), lastImageInt(i,8)];
                                LIFFTPCounter = LIFFTPCounter + 1;
                            end
                        end
                    else
                        if(lastImageInt(i,1) <= size(im1,2) && lastImageInt(i,1) >= 1 && lastImageInt(i,2) <= size(im1,1) && lastImageInt(i,2) >= 1)
                            firstFoundTrackingPoints(FFTPCounter, :) = [round(lastImageInt(i,1)), round(lastImageInt(i,2)), imageNumber, lastImageInt(i,4), guessedParticleTrack, lastImageInt(i,7), lastImageInt(i,8)];
                            FFTPCounter = FFTPCounter + 1;
                            lastImageFFTP(LIFFTPCounter, :) = [round(lastImageInt(i,1)), round(lastImageInt(i,2)), imageNumber+1, lastImageInt(i,4), guessedParticleTrack, lastImageInt(i,7), lastImageInt(i,8)];
                            LIFFTPCounter = LIFFTPCounter + 1;
                        end
                    end
                end
            end
            coord2 = coord2(coord2(:,3) == 0,:);
            %We are adding all the new particles that could not be
            %connected to particles in the previous image.
            for j=1:size(coord2,1)
                lastImageFFTP(LIFFTPCounter, :) = [coord2(j,1), coord2(j,2), imageNumber+1, uniqueID, keepParticleTrack, Uint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Uint(coord2(j,2), coord2(j,1))*previousVelocityContribution, Vint(coord2(j,2), coord2(j,1))*(1-previousVelocityContribution) + Vint(coord2(j,2), coord2(j,1))*previousVelocityContribution];
                LIFFTPCounter = LIFFTPCounter + 1;
                uniqueID = uniqueID + 1;
            end
        end
        %Removing excess zero rows that were not filled.
        lastImageFFTP(lastImageFFTP(:,4)==0,:) = [];
    end
    %Removing excess zero rows that were not filled.
    firstFoundTrackingPoints(firstFoundTrackingPoints(:,4)==0,:) = [];
    IDsortedTracks = sortrows(firstFoundTrackingPoints,4);
    %Removing too many guessed values.
    if(particleSkippedFrames > 0)
        findDifferenceIndexID = find(diff(IDsortedTracks(:,4)) ~= 0);
        trajectoryLengthsInArray = diff(findDifferenceIndexID);
        guessCounter = 0;
        keepCounter = 0;
        firstGuessPosition = 0;
        %If there are more particles with a given ID than
        %particleSkippedFrames then we will check it to see if it is needed
        %to remove the guessed particle tracks.
        if(~isempty(findDifferenceIndexID))
            if(findDifferenceIndexID(1) > particleSkippedFrames)
                sliceOfID = IDsortedTracks(1:findDifferenceIndexID(1),:);
                for j=1:size(sliceOfID)
                    if(sliceOfID(j,5) == 1)
                        if(firstGuessPosition == 0)
                            firstGuessPosition = j;
                        end
                        guessCounter = guessCounter + 1;
                        if(j == size(sliceOfID))
                            if(guessCounter > particleSkippedFrames)
                                IDsortedTracks(firstGuessPosition:j,:) = 0;
                                keepCounter = 0;
                            end
                        end
                    else
                        if(guessCounter > particleSkippedFrames)
                            IDsortedTracks(firstGuessPosition:j,:) = 0;
                            keepCounter = 0;
                        else
                            keepCounter = keepCounter + 1;
                        end
                        firstGuessPosition = 0;
                    end
                end
            end
            if(size(IDsortedTracks,1) - findDifferenceIndexID(end) > particleSkippedFrames)
                sliceOfID = IDsortedTracks(findDifferenceIndexID(end)+1:size(IDsortedTracks,1),:);
                for j=1:size(sliceOfID)
                    if(sliceOfID(j,5) == 1)
                        if(firstGuessPosition == 0)
                            firstGuessPosition = j;
                        end
                        guessCounter = guessCounter + 1;
                        if(j == size(sliceOfID))
                            if(guessCounter > particleSkippedFrames)
                                IDsortedTracks(firstGuessPosition:j,:) = 0;
                                keepCounter = 0;
                            end
                        end
                    else
                        if(guessCounter > particleSkippedFrames)
                            IDsortedTracks(firstGuessPosition:j,:) = 0;
                            keepCounter = 0;
                        else
                            keepCounter = keepCounter + 1;
                        end
                        firstGuessPosition = 0;
                    end
                end
            end
            for i=1:size(trajectoryLengthsInArray,1)-1
                if(trajectoryLengthsInArray(i) > particleSkippedFrames)
                    sliceOfID = IDsortedTracks(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1),:);
                    for j=1:size(sliceOfID)
                        if(sliceOfID(j,5) == 1)
                            if(firstGuessPosition == 0)
                                firstGuessPosition = j;
                            end
                            guessCounter = guessCounter + 1;
                            if(j == size(sliceOfID))
                                if(guessCounter > particleSkippedFrames)
                                    IDsortedTracks(firstGuessPosition:j,:) = 0;
                                    keepCounter = 0;
                                end
                            end
                        else
                            if(guessCounter > particleSkippedFrames)
                                IDsortedTracks(firstGuessPosition:j,:) = 0;
                                keepCounter = 0;
                            else
                                keepCounter = keepCounter + 1;
                            end
                            firstGuessPosition = 0;
                        end
                    end
                end
            end
            IDsortedTracks(IDsortedTracks(:,4)==0,:) = [];
        end
    end
    
    %Extracting the start and end-point of every trajectory.
    trajectoryEndPoints = IDsortedTracks;
    findDifferenceIndexID = find(diff(trajectoryEndPoints(:,4)) ~= 0);
    if(~isempty(findDifferenceIndexID))
        trajectoryLengthsInArray = diff(findDifferenceIndexID);
        if(findDifferenceIndexID(1) >= 2)
            trajectoryEndPoints(2:findDifferenceIndexID(1)-1,:) = 0;
        else
            trajectoryEndPoints(1:findDifferenceIndexID(1),:) = 0;
        end
        if(size(trajectoryEndPoints,1) - findDifferenceIndexID(end) >= 2)
            trajectoryEndPoints(findDifferenceIndexID(end)+2:size(trajectoryEndPoints,1)-1,:) = 0;
        else
            trajectoryEndPoints(findDifferenceIndexID(end)+1:size(trajectoryEndPoints,1),:) = 0;
        end
        for i=1:size(trajectoryLengthsInArray,1)
            if(trajectoryLengthsInArray(i) >= 2)
                trajectoryEndPoints(findDifferenceIndexID(i)+2:findDifferenceIndexID(i+1)-1,:) = 0;
            else
                trajectoryEndPoints(findDifferenceIndexID(i)+1:findDifferenceIndexID(i+1),:) = 0;
            end
        end
        trajectoryEndPoints(trajectoryEndPoints(:,4)==0,:) = [];

        %Connect trajectories that have start-coordinate within the given 
        %trackingSearchRadius of a end-coordinate of another trajectory, and that 
        %are in chronological order
        firstFoundTrackingPoints = IDsortedTracks;
        for i=1:size(trajectoryEndPoints,1)-1
            for j=i+1:size(trajectoryEndPoints,1)
                if((trajectoryEndPoints(i,1)-trajectoryEndPoints(j,1))^2 + (trajectoryEndPoints(i,2)-trajectoryEndPoints(j,2))^2 <= trackingSearchRadiusSquared)
                    if((trajectoryEndPoints(i,3) + 1) == trajectoryEndPoints(j,3))
                        firstFoundTrackingPoints(firstFoundTrackingPoints(:,4) == trajectoryEndPoints(j,4),4) = trajectoryEndPoints(i,4);
                    end
                end
            end
        end
    end
    
    try
        close(waitBar)
    catch
        disp("No loading bar to close. Exiting PTV")
    end
    
    %return
    if nargout == 1
        varargout{1} = firstFoundTrackingPoints;
    elseif nargout == 2
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = piv;
    elseif nargout == 3
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = piv;
        varargout{3} = Uint;
    elseif nargout == 4
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = piv;
        varargout{3} = Uint;
        varargout{4} = Vint;
    elseif nargout == 5
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = piv;
        varargout{3} = Uint;
        varargout{4} = Vint;
        varargout{5} = ensembleUint;
    elseif nargout == 6
        varargout{1} = firstFoundTrackingPoints;
        varargout{2} = piv;
        varargout{3} = Uint;
        varargout{4} = Vint;
        varargout{5} = ensembleUint;
        varargout{6} = ensembleVint;
    end
end