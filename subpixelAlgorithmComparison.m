p = genpath('HydrolabPIV\src');
addpath(p);
javaaddpath(['HydrolabPIV\src' filesep 'measures']);
javaaddpath(['HydrolabPIV\src' filesep 'interp']);
javaaddpath(['HydrolabPIV\src' filesep 'subpixel']);
p = genpath('HydrolabPIV');
addpath(p);
javaaddpath(['HydrolabPIV' filesep 'images']);

%Testing speed and accuracy of subpixel functions.
im1 = im2double(imread('imA.png'));
sizeOfWindow = [32,64,128,256,500];
dxTrue = [5;4];
t1 = []; t2 = []; t3 = []; t4 = []; t5 = []; t6 = [];
e1 = []; e2 = []; e3 = []; e4 = []; e5 = []; e6 = [];
t1Mean = []; t2Mean = []; t3Mean = []; t4Mean = []; t5Mean = []; t6Mean = [];
e1Mean = []; e2Mean = []; e3Mean = []; e4Mean = []; e5Mean = []; e6Mean = [];
numberOfObservations = 1000;
for i=sizeOfWindow
    for j=1:numberOfObservations
        idx = (1:i);
        A = im1(idx,idx);
        B = im1(idx+dxTrue(2),idx+dxTrue(1));

        % Calculate cross-correlation
        xm = i/2;
        idy = (-xm:xm) + i; 
        ncc = standardcc(A,[],B,[],idy,idy);

        % Find displacement dx
        tic
        [x0,delta,out] = subpixel3x2(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t1 = [t1,toc];
        e1 = [e1,mean(dxTrue-dx)];
        tic
        [x0,delta,out] = subpixel3x3(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t2 = [t2,toc];
        e2 = [e2,mean(dxTrue-dx)];
        tic
        [x0,delta,out] = subpixel3x3lm(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t3 = [t3,toc];
        e3 = [e3,mean(dxTrue-dx)];
        tic
        [x0,delta,out] = subpixel3x3ls(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t4 = [t4,toc];
        e4 = [e4,mean(dxTrue-dx)];
        tic
        [x0,delta,out] = subpixel5x5lm(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t5 = [t5,toc];
        e5 = [e5,mean(dxTrue-dx)];
        tic
        [x0,delta,out] = subpixel5x5ls(ncc(idy,idy));
        dx = x0-delta-xm-1;
        t6 = [t6,toc];
        e6 = [e6,mean(dxTrue-dx)];
    end
    
    t1Mean = [t1Mean, mean(t1)];
    e1Mean = [e1Mean, mean(e1)];
    t1 = [];
    e1 = [];
    t2Mean = [t2Mean, mean(t2)];
    e2Mean = [e2Mean, mean(e2)];
    t2 = [];
    e2 = [];
    t3Mean = [t3Mean, mean(t3)];
    e3Mean = [e3Mean, mean(e3)];
    t3 = [];
    e3 = [];
    t4Mean = [t4Mean, mean(t4)];
    e4Mean = [e4Mean, mean(e4)];
    t4 = [];
    e4 = [];
    t5Mean = [t5Mean, mean(t5)];
    e5Mean = [e5Mean, mean(e5)];
    t5 = [];
    e5 = [];
    t6Mean = [t6Mean, mean(t6)];
    e6Mean = [e6Mean, mean(e6)];
    t6 = [];
    e6 = [];
end
fontSize = 30;
figure
plot(t1Mean,e1Mean,'O',t2Mean,e2Mean,'O',t3Mean,e3Mean,'O',t4Mean,e4Mean,'O',t5Mean,e5Mean,'O',t6Mean,e6Mean,'O', 'LineWidth', 12)
legend('3x2', '3x3', '3x3lm', '3x3ls', '5x5lm','5x5ls')
titleText = ['Comparison of subpixel algorithms: N=', num2str(numberOfObservations)];
title(titleText)
xlabel('Time taken in seconds (s)') % x-axis label
ylabel('Error') % y-axis label
set(gca, 'FontSize', fontSize)