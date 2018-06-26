function [varargout]=drawCircles(image,particlePositions,particleDiameter,blur,sigma)
    %Minimum 3 inputs and maximum 5 inputs
    narginchk(3,5)
    if ~exist('blur', 'var')
        blur = 'none';
        radius = particleDiameter/2;
    else
        radius = particleDiameter/2 - 1;
        if(radius < 0.5)
            radius = 0.5;
        end
    end
    if ~exist('sigma', 'var')
        sigma = 0.5;
    end
    
    imageZeros = zeros(size(image));
    particlePositions = [particlePositions, ones(size(particlePositions,1),1)*radius];
    imageWithParticles = insertShape(imageZeros, 'FilledCircle', particlePositions, 'color', 'white', 'opacity', 1);
    
    if(strcmp(blur,'none'))
        if(size(imageWithParticles,3) == 3)
            varargout{1} = im2double(rgb2gray(imageWithParticles));
        else
            varargout{1} = im2double(imageWithParticles);
        end
    else
        if(size(imageWithParticles,3) == 3)
            blurredImage = rgb2gray(imgaussfilt(imageWithParticles, sigma));
            varargout{1} = im2double(blurredImage);
        else
            blurredImage = imgaussfilt(imageWithParticles, sigma);
            varargout{1} = im2double(blurredImage);
        end
    end
end