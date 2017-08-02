function [ result ] = DrawBox( frm, Height, Width, x, y, Rx, Ry)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    result = frm;
    
    ist = max(1, floor(x-Rx));
    ied = min(Height, ceil(x+Rx));
    jst = max(1, floor(y-Ry));
    jed = min(Width, ceil(y+Ry));

    for i=ist:ied
        result(i, jst+0, :) = [1, 0, 0];
        if (jst+1 <= Width)
            result(i, jst+1, :) = [1, 0, 0];
        end
        if (jst+2 <= Width)
            result(i, jst+2, :) = [1, 0, 0];
        end
        result(i, jed-0, :) = [1, 0, 0];
        if (jed-1 >= 1)
            result(i, jed-1, :) = [1, 0, 0];
        end
        if (jed-2 >= 1)
            result(i, jed-2, :) = [1, 0, 0];
        end
    end
    for j=jst:jed
        result(ist+0, j, :) = [1, 0, 0];
        if (ist+1 <= Height)
            result(ist+1, j, :) = [1, 0, 0];
        end
        if (ist+2 <= Height)
            result(ist+2, j, :) = [1, 0, 0];
        end
        result(ied-0, j, :) = [1, 0, 0];
        if (ied-1 >= 1)
            result(ied-1, j, :) = [1, 0, 0];
        end
        if (ied-2 >= 1)
            result(ied-2, j, :) = [1, 0, 0];
        end
    end
end

