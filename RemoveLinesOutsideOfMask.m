function [sortedFFTP]=RemoveLinesOutsideOfMask(sortedFFTPOriginal, mask, N)
    sortedFFTP = sortedFFTPOriginal;
    if(mask == 0)
        try
            maskFile = uigetfile(pwd, 'Select Mask', '*.*');
        catch
            disp('No mask selected. Exiting RemoveLinessOutsideOfMask function.')
            return
        end
        if(ischar(maskFile) == 0)
            disp('No mask selected. Exiting RemoveLinessOutsideOfMask function.')
            return
        else
            mask = bwdist(imcomplement(imread(maskFile)));
        end
    end
    if(N <= 1)
        N = 2;
    end
    for i=1:size(sortedFFTPOriginal,1)-1
        if(sortedFFTPOriginal(i,4) == sortedFFTPOriginal(i+1,4))
            if(sortedFFTPOriginal(i,1) ~= sortedFFTPOriginal(i+1,1))
                a = (sortedFFTPOriginal(i+1,2) - sortedFFTPOriginal(i,2))/(sortedFFTPOriginal(i+1,1) - sortedFFTPOriginal(i,1));
                b = sortedFFTPOriginal(i,2) - sortedFFTPOriginal(i,1)*a;
                deltaX = (sortedFFTPOriginal(i+1,1) - sortedFFTPOriginal(i,1))/N;
                for j=1:N
                    nextX = round(sortedFFTPOriginal(i,1)+deltaX*j);
                    nextY = round(nextX*a + b);
                    if(nextX < 1)
                        nextX = 1;
                    end
                    if(nextY < 1)
                        nextY = 1;
                    end
                    if(mask(nextY,nextX) == 0)
                        sortedFFTP(sortedFFTP(:,4) == sortedFFTP(i,4),:) = 0;
                        break
                    end
                end
            end
        end
    end
    sortedFFTP(sortedFFTP(:,4) == 0,:) = [];
end