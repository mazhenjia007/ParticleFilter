function [ lnp ] = CalLikelihood_Subpixel( frm, Rad, is, js )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [Height, Width] = size(frm); 
    
    n0A = sum(sum(frm==uint8(0)));
    n1A = sum(sum(frm==uint8(1)));
    
    n0D = 0; % for correct
    n1D = 0; % for incorrect
    
    ist = floor(is - Rad + 0.5); ist = max(ist, 1);
    ied = ceil(is + Rad + 0.5); ied = min(ied, Height);
    jst = floor(js - Rad + 0.5); jst = max(jst, 1);
    jed = ceil(js + Rad + 0.5); jed = min(jed, Width);

    for i=ist:ied
        for j=jst:jed
%     for i=1:Height
%         for j=1:Width
            dlt = (is-(i-0.5))^2 + (js-(j-0.5))^2 - Rad^2;
            if dlt < 0  % in
                if frm(i, j) == uint8(0)
                    n0D = n0D + 1;
                else
                    n1D = n1D + 1;
                end
            end
        end
    end
    
    n0B = n0A - n0D;
    n1B = n1A - n1D;

    p_white = n1B / (n0B + n1B);
    % p_white = min(p_white+0.3, 0.95);
    p_white = max(p_white-0.3, 0.03);
    lnp = (n0D)*log(p_white) + (n1D)*log(1-p_white);
    lnp = lnp / (n0D + n1D);
    % lnp = lnp / abs(log((n0A)/(n0A+n1A)));
end

