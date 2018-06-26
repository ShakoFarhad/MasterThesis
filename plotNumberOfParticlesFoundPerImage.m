function [] = plotNumberOfParticlesFoundPerImage(images, threshold, particleLength, filter, noiseLength)

    numberOfParticles = [];
    frames = 1:size(images,2);
    for imageNumber=1:size(images,2)
        tempImage1 = imread(images{imageNumber});
        if(size(tempImage1,3) == 3)
            im1 = im2double(rgb2gray(tempImage1));
        else
            im1 = im2double(tempImage1(:,:,1));
        end

        coordinates = findParticles(im1, threshold, particleLength, filter, noiseLength);
        numberOfParticles = [numberOfParticles, size(coordinates,1)];
    end
    p = polyfit(frames,numberOfParticles,1);
    fontSize = 30;
    figure
    plot(frames,numberOfParticles,'b', [frames(1), frames(end)], [p(1)*frames(1) + p(2),p(1)*frames(end) + p(2)], 'r-', 'LineWidth', 1.5);
    legend('Particle Data', ['Linear Fit: y = ', num2str(p(1)), 'x + ', num2str(p(2))])
    titleText = ['Number of Particles per Frame. Threshold: ', num2str(threshold), '.'];
    title(titleText)
    xlabel('Frame Number') % x-axis label
    ylabel('Number of Particles') % y-axis label
    set(gca, 'FontSize', fontSize)
end

