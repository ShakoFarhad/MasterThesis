p = genpath('HydrolabPIV\src');
addpath(p);
javaaddpath(['HydrolabPIV\src' filesep 'measures']);
javaaddpath(['HydrolabPIV\src' filesep 'interp']);

AmountOfParticles = 300;
particleSpawnRate = 2;
BGX = 512;
BGY = 512;
emitXArray = [30, 300];%[1,70];
emitYArray = [480, 512];%[1,70];
%VX = zeros([BGX, BGY]);
%VY = zeros([BGX, BGY]);
localXVelocity = 1.5;
localYVelocity = 1.5;

mu = [0 0]; Sigma = [2 1; 1 2];
x1 = linspace(-2,2,BGX); x2 = linspace(-2,2,BGY);
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(BGX+BGY)/8;
VXX = reshape(F,length(x2),length(x1));
VYY = -reshape(F,length(x2),length(x1));
VXX = rot90(imgaussfilt(VXX,700));

particleDiameter = [9, 15];
particleLumi0 = [0.5,1.0];

myMask = imcomplement(gray);
subWindow = 30;
searchArea = 10;
overlap = 0.5;
multipass = [];
imageFolder = 'G:\2016-08-11_10_10.44.34_1GB'; %Hjemme
addpath(imageFolder);
image = dir(imageFolder);
fieldNames = fieldnames(image);
images = rmfield(image, fieldNames(2:end));
images = struct2cell(images);
images(:,1:3) = [];
masks = {'G:\Master_Thesis_2017-2018\Code\Images\PIVMaskFolder\InvertedDistanceImageSelfMade.png'};
mask1 = imread(masks{1});

im1PIV = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
im2PIV = zeros(size(imread(images{1}),1),size(imread(images{1}),2),floor(size(images,2)/2));
n = 1;
for imageNumber=1:2:size(images,2)-1
    tempImage1 = imread(images{imageNumber});
    tempImage2 = imread(images{imageNumber + 1});
    im1PIV(:,:,n) = im2double(tempImage1(:,:,1));
    im2PIV(:,:,n) = im2double(tempImage2(:,:,1));
    n = n + 1;
end
opt = setpivopt('range',[-searchArea searchArea -searchArea searchArea],'subwindow',subWindow,subWindow,overlap,'ensemble',@nanmedian);
piv = normalpass(multipass,im1PIV,mask1,im2PIV,mask1,opt);
[VX,VY] = replaceoutliers(piv);

x = linspace(0, size(VX,1), size(VX,1));
y = linspace(0, size(VX,2), size(VX,2));
xq = linspace(0, size(VX,1), size(im1PIV,1));
yq = linspace(0, size(VX,2), size(im1PIV,2));
[X, Y] = meshgrid(x,y);
[Xq, Yq] = meshgrid(xq,yq);
Uint = interp2(X,Y,VX,Xq,Yq, 'cubic');

x = linspace(0, size(VY,1), size(VY,1));
y = linspace(0, size(VY,2), size(VY,2));
xq = linspace(0, size(VY,1), size(im1PIV,1));
yq = linspace(0, size(VY,2), size(im1PIV,2));
[X, Y] = meshgrid(x,y);
[Xq, Yq] = meshgrid(xq,yq);
Vint = interp2(X,Y,VY,Xq,Yq, 'cubic');

Uint = flipud(Uint);
Vint = flipud(Vint);

blurU = imgaussfilt(Uint, 20);
blurV = imgaussfilt(Vint, 20);

figure
surf(xq, yq, blurV)
xlabel('x')
ylabel('y')

figure
load wind
pcolor(xq,yq,Vint); shading interp
hold on;
d = quiver(xq, yq, Uint, Vint, 3, 'w');
hold off
title(['Velocity field. Subwindow: ', num2str(subWindow), '. Search area: ', num2str(searchArea)])
xlabel('Pixels') % x-axis label
ylabel('Pixels') % y-axis label
colormap winter
set(gca, 'YDir', 'reverse')
c = colorbar;

particleEmitter = ParticleEmitter(AmountOfParticles, particleSpawnRate, BGX, BGY, emitXArray, emitYArray, VXX, VYY, localXVelocity, localYVelocity);
particleEmitter.initializeParticles(particleDiameter, particleLumi0);
PA0 = particleEmitter.ParticleArray;
amountOfUpdates = 120;
saveFolder = uigetdir('C:\', 'Select Synthetic Image Save Folder');
for k=1:round(AmountOfParticles/particleSpawnRate+amountOfUpdates)
    if(mod(k,10))
        disp([num2str(k/(AmountOfParticles/particleSpawnRate+amountOfUpdates)*100), '%'])
    end
    [window, PA] = update(particleEmitter);
    figure('visible', 'off');
    imagesc(window)
    imageName = [saveFolder, sprintf('MySavedPlotNoTracks%d.png', k)];
    imwrite(window,imageName);
    
%     hold on
%     i = 1;
%     while(i <= AmountOfParticles)
%        plot(cell2mat(PA(i).YTrack), cell2mat(PA(i).XTrack), 'LineWidth', 2)
%        i = i + 1;
%     end
%     axis([1 BGX 1 BGY])
%     colormap('gray')
%     ax = gca;
%     ax.Visible = 'off';
%     outerpos = ax.OuterPosition;
%     left = outerpos(1);
%     bottom = outerpos(2);
%     ax_width = outerpos(3);
%     ax_height = outerpos(4);
%     ax.Position = [left bottom ax_width ax_height];
%     imageName = sprintf('SyntheticParticleData\\WithTracks\\MySavedPlotTracks%d.png', k);
%     f=getframe; imwrite(f.cdata,imageName);
end

%Det under er bare tull for aa teste stuff
% imageFolder = char('G:\2016-08-11_10_10.44.34_1GB');
% addpath(imageFolder);
% image = dir(imageFolder);
% fieldNames = fieldnames(image);
% images = rmfield(image, fieldNames(2:end));
% images = struct2cell(images);
% images(:,1:2) = [];
% tempImage = imread(images{1});
% image = im2double(tempImage(:,:,1));
% particleDiameter = 4;
% coordinates = findParticles(image, 5.5, particleDiameter);
% figure
% imshow(image)
% figure
% imshow(drawCircles(image, coordinates, particleDiameter, 'gaussian', 0.5))