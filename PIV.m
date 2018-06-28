function [varargout] = PIV(im1, im2, mask1, mask2, overlap, subWindow, searchArea, replaceOutliers, distortedPass, multipass, subpixel, particleLength)

%     p = genpath('HydrolabPIV\src');
%     addpath(p);
%     javaaddpath(['HydrolabPIV\src' filesep 'measures']);
%     javaaddpath(['HydrolabPIV\src' filesep 'interp']);
    
    %Minimum 2 inputs and maximum 12 inputs
    narginchk(2,12)
    
    if ~exist('mask1', 'var')
        mask1 = [];
    end
    if ~exist('mask2', 'var')
        mask2 = [];
    end
    if ~exist('overlap', 'var')
        overlap = 0.5;
    end
    if ~exist('subWindow', 'var')
        imageSize = max(size(im1));
        subWindow = floor(imageSize/10);
    end
    if ~exist('searchArea', 'var')
        if(size(subWindow,2) == 1)
            searchArea = floor(subWindow(1)/3);
        else
            searchArea = floor(min(subWindow)/3);
        end
    end
    if ~exist('replaceOutliers', 'var')
        replaceOutliers = 2;
    end
    if ~exist('distortedPass', 'var')
        distortedPass = 0;
    end
    if ~exist('multipass', 'var')
        multipass = [];
    end
    if ~exist('subpixel', 'var')
        subpixel = @subpixelnone;
    end
    if ~exist('particleLength', 'var')
        particleLength = 0;
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
    
    
    if(distortedPass)
        if(particleLength > 4)
            opt = setpivopt('iminterp',@iminterp2bsplinej,[],'subpixel', subpixel);
        elseif(particleLength <= 4)
            opt = setpivopt('iminterp',@iminterp2lanczosj,3,'subpixel', subpixel);
        elseif(particleLength == 0)
            opt = setpivopt('iminterp',@iminterp2matlab,'linear','subpixel', subpixel);
        end
        piv = distortedpass(multipass,im1,mask1,im2,mask2,opt);
    else
        opt = setpivopt('range',[-x1 x1 -y1 y1],'subwindow',x,y,overlap,'subpixel', subpixel);
        piv = normalpass(multipass,im1,mask1,im2,mask2,opt);
    end
    
    if replaceOutliers == 1
        if isempty(mask1) || isempty(mask2)
            [U,V,x,y] = replaceoutliers(piv);
            piv.U = U;
            piv.V = V;
        else
            [U,V,x,y] = replaceoutliers(piv,mask1&mask2);
            piv.U = U;
            piv.V = V;
        end
    elseif replaceOutliers == 2
        [U,V,x,y] = replaceoutliers(piv);
        piv.U = U;
        piv.V = V;
    else
        U = piv.U;
        V = piv.V;
    end
    
    %return
    if nargout == 1
        varargout{1} = piv;
    elseif nargout == 2
        varargout{1} = piv;
        varargout{2} = U;
    elseif nargout == 3
        varargout{1} = piv;
        varargout{2} = U;
        varargout{3} = V;
    elseif nargout == 4
        varargout{1} = piv;
        varargout{2} = U;
        varargout{3} = V;
        varargout{4} = x;
    elseif nargout == 4
        varargout{1} = piv;
        varargout{2} = U;
        varargout{3} = V;
        varargout{4} = x;
        varargout{5} = y;
    end
end