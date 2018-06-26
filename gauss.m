function [G] = gauss(rows, columns, sigma)
    [x, y] = meshgrid((round(-rows/2)+1):(round(rows/2)-1), (round(-columns/2)+1):(round(columns/2)-1));
    %If sigma is not 10 times smaller than the mesh size then
    %the sum of G will be less than 1 because significant values
    %gets cut off. Therefore we need to normalize by dividing by
    %the sum of G, and this is significantly slower for big meshes.
    if(10*sigma <= min(rows,columns))
        G = 1/(2*pi*sigma^2)*exp(-x.^2/(2*sigma^2) - y.^2/(2*sigma^2));
    else
        G = exp(-x.^2/(2*sigma^2) - y.^2/(2*sigma^2));
        G./sum(G(:));
    end
end