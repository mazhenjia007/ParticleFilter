function [ lnp ] = CalLikelihood( frm, ptn, i, j, Alpha0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [hF, wF] = size(frm);
    [hP, wP] = size(ptn);
    
    istf = i; iedf = i+1;
    jstf = j; jedf = j+1;
    istp = hP/2; iedp = hP/2+1;
    jstp = wP/2; jedp = wP/2+1;
    
    while (istf>0 && istp>0)
        istf = istf - 1;
        istp = istp - 1;
    end
    istf = istf + 1; istp = istp + 1;
    
    while (iedf<=hF && iedp<=hP)
        iedf = iedf + 1;
        iedp = iedp + 1;
    end
    iedf = iedf - 1; iedp = iedp - 1;
    
    while (jstf>0 && jstp>0)
        jstf = jstf - 1;
        jstp = jstp - 1;
    end
    jstf = jstf + 1; jstp = jstp + 1;
    
    while (jedf<=wF && jedp<=wP)
        jedf = jedf + 1;
        jedp = jedp + 1;
    end
    jedf = jedf - 1; jedp = jedp - 1;

    n0F_all = sum(sum(frm==0));
    n1F_all = sum(sum(frm==1));
    
%     n0P = sum(ptn(istp:iedp, jstp:jedp)==0);
%     n1P = sum(ptn(istp:iedp, jstp:jedp)==1);
    
    n0F = sum(sum(frm(istf:iedf, jstf:jedf)==0));
    n1F = sum(sum(frm(istf:iedf, jstf:jedf)==1));
    
    n0B = n0F_all - n0F;
    n1B = n1F_all - n1F;
    
    p_noise = n0B / (n0B + n1B);
    
    dif = abs(frm(istf:iedf, jstf:jedf)-ptn(istp:iedp, jstp:jedp));
    
    n0D = sum(sum(dif==0)); % for correct pix
    n1D = sum(sum(dif==1)); % for incorrect pix
    
    lnp = (n0D+Alpha0)*log(1-p_noise) + (n1D+Alpha0)*log(p_noise);
    lnp = lnp / numel(dif);
end

