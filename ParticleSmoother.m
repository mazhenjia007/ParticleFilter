load('img400b.mat');

IMHEIGHT = 300;
IMWIDTH = 400;

eps = 1e-7;

[Height, Width, nFrame] = size(vs);

% Set pattern
hPattern = 6; wPattern = 6;
fPattern = [1, 1, 0, 0, 1, 1;       % center (3.5, 3.5)
            1, 0, 0, 0, 0, 1;
            0, 0, 0, 0, 0, 0;
            0 ,0 ,0 ,0 ,0, 0;
            1, 0, 0, 0, 0, 1;
            1, 1, 0, 0, 1, 1];
fPattern = uint8(fPattern);
Rad = 3;

% #particles
M = 7000;

% [x1, x2, v1, v2] : x1 v1 in Height direction, x2 v2 in Width direction

sigmax = 0.05;
sigmav = 0.1;

p = zeros(M, 4);
w = zeros(M, 1);

p(:, 1) = rand(M, 1) * Height;
p(:, 2) = rand(M, 1) * Width;

ps = zeros(M, 4, nFrame);
ws = zeros(M, nFrame);

% particles = SetInitParticles(M, Height, Width);

% figId = figure(1);

hwb = waitbar(0, 'Initializing...');

tag = 0;

for iFrame = 1:nFrame
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Forward: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));

    % Dynamic   
    p = ParticleDynamic(p, M, Height, Width, sigmax, sigmav);
        
    % Weighting ( Likelihood )
    for m=1:M
        w(m) = CalLikelihood_Subpixel(vs(:, :, iFrame), ...
            Rad, p(m, 1), p(m, 2));
        w(m) = exp(w(m));
    end
    
    w = w / sum(w);
    
    ps(:, :, iFrame) = p;
    ws(:, iFrame) = w;
    
    % Resampling
    plabel = randsample(1:M, M, true, w);
    p = p(plabel, :);
end

% wT = zeros(M, 1);
wTs = zeros(M, nFrame);

wT = w;
wT = wT / sum(wT);

wTs(:, nFrame) = wT;

tag = 0;
for iFrame=nFrame-1:-1:1
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Backward: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));
        
    dorm = zeros(M, 1);
    for k=1:M
        pkx = ps(k, 1:2, iFrame) + ps(k, 3:4, iFrame);
        pkv = ps(k, 3:4, iFrame);
        
        pjx = ps(:, 1:2, iFrame+1);
        pjv = ps(:, 3:4, iFrame+1);
        
        dx = bsxfun(@minus, pjx, pkx) / sigmax;
        dv = bsxfun(@minus, pjv, pkv) / sigmav;
        
        dlt = (- dx(:, 1).^2 - dx(:, 2).^2 - dv(:, 1).^2 - dv(:, 2).^2) / 2;
        
        pdf = bsxfun(@times, ws(:, iFrame), exp(dlt));
        dorm(:, 1) = dorm(:, 1) + pdf;
    end
    
    COEj = zeros(M, 1);
    for j=1:M
        if (dorm(j, 1) < eps)
            COEj(j, 1) = 0;
        else
            if (wTs(j, iFrame+1) < eps)
                % COEj(j, 1) = 1e-10;
                COEj(j, 1) = 0;
            else
                COEj(j, 1) = wTs(j, iFrame+1) / dorm(j, 1);
            end
        end
    end
    
%     for j=1:M
%         if (wTs(j, iFrame+1) < eps)
%             COEj(j, 1) = 0;
%         else
%             COEj(j, 1) = wTs(j, iFrame+1) / dorm(j, 1);
%         end
%     end
    
    wT = zeros(M, 1);
    for i=1:M
        pix = ps(i, 1:2, iFrame) + ps(i, 3:4, iFrame);
        piv = ps(i, 3:4, iFrame);
        
        pjx = ps(:, 1:2, iFrame+1);
        pjv = ps(:, 3:4, iFrame+1);
        
        dx = bsxfun(@minus, pjx, pix) / sigmax;
        dv = bsxfun(@minus, pjv, piv) / sigmav;
        
        dlt = (- dx(:, 1).^2 - dx(:, 2).^2 - dv(:, 1).^2 - dv(:, 2).^2) / 2;
        
        wT(i) = ws(i, iFrame) * sum(exp(dlt).*COEj(:, 1));
        
        if (wT(i) < eps || isnan(wT(i)))
            wT(i) = 1e-10;
        end
    end
    
    wT = wT / sum(wT);
    
    wTs(:, iFrame) = wT;
end

impFrame = zeros(IMHEIGHT, IMWIDTH, nFrame);

lastx1 = -1;
lasty1 = -1;
lastx2 = -1;
lasty2 = -1;
tjd1Frame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
tjd2Frame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
lastFrame1 = zeros(IMHEIGHT, IMWIDTH, 3);
lastFrame2 = zeros(IMHEIGHT, IMWIDTH, 3);

mus_f = zeros(4, nFrame);
sigs_f = zeros(4, nFrame);
mus_fb = zeros(4, nFrame);
sigs_fb = zeros(4, nFrame);

mnpFrame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);

tag = 0;
for iFrame=1:nFrame
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Visual: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));

    imp = zeros(IMHEIGHT, IMWIDTH);
    for m=1:M
        imp(ceil(ps(m, 1, iFrame)/Height*IMHEIGHT), ceil(ps(m, 2, iFrame)/Width*IMWIDTH)) = 1;
    end
    impFrame(:, :, iFrame) = imp;
    
    mx1 = sum(ps(:, 1, iFrame).*ws(:, iFrame)) / sum(ws(:, iFrame)) / Height * IMHEIGHT;
    my1 = sum(ps(:, 2, iFrame).*ws(:, iFrame)) / sum(ws(:, iFrame)) / Width * IMWIDTH;

    if (lastx1 ~= -1)
        dx = mx1 - lastx1;
        dy = my1 - lasty1;
        dt = min(floor(abs(dx/Height*IMHEIGHT)), floor(abs(dy/Width*IMWIDTH)));
        dx = dx / dt;
        dy = dy / dt;
        for t=1:dt
            lastFrame1 = DrawBox(lastFrame1, IMHEIGHT, IMWIDTH, lastx1+dx*t, lasty1+dy*t, 0, 0);
        end
    end
    lastFrame1 = DrawBox(lastFrame1, IMHEIGHT, IMWIDTH, mx1, my1, 0, 0);
    tjd1Frame(:, :, :, iFrame) = lastFrame1;
    lastx1 = mx1;
    lasty1 = my1;
    
    mx2 = sum(ps(:, 1, iFrame).*wTs(:, iFrame)) / sum(wTs(:, iFrame)) / Height * IMHEIGHT;
    my2 = sum(ps(:, 2, iFrame).*wTs(:, iFrame)) / sum(wTs(:, iFrame)) / Width * IMWIDTH;

    if (lastx2 ~= -1)
        dx = mx2 - lastx2;
        dy = my2 - lasty2;
        dt = min(floor(abs(dx/Height*IMHEIGHT)), floor(abs(dy/Width*IMWIDTH)));
        dx = dx / dt;
        dy = dy / dt;
        for t=1:dt
            lastFrame2 = DrawBox(lastFrame2, IMHEIGHT, IMWIDTH, lastx2+dx*t, lasty2+dy*t, 0, 0);
        end
    end
    lastFrame2 = DrawBox(lastFrame2, IMHEIGHT, IMWIDTH, mx2, my2, 0, 0);
    tjd2Frame(:, :, :, iFrame) = lastFrame2;
    lastx2 = mx2;
    lasty2 = my2;
    
    
    mnp = ind2rgb(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), [0, 0, 0; 1, 1, 1]);
    mnp = DrawBox(mnp, IMHEIGHT, IMWIDTH, mx2, my2, Rad/Height*IMHEIGHT, Rad/Width*IMWIDTH);
    mnpFrame(:, :, :, iFrame) = mnp;
    

    mu_f = zeros(4, 1);
    sig_f = zeros(4, 1);
    mu_fb = zeros(4, 1);
    sig_fb = zeros(4, 1);
    for m=1:M
        mu_f = mu_f + ws(m, iFrame) * ps(m, :, iFrame)';
        sig_f = sig_f + ws(m, iFrame) * (ps(m, :, iFrame).^2)';
        mu_fb = mu_fb + wTs(m, iFrame) * ps(m, :, iFrame)';
        sig_fb = sig_fb + wTs(m, iFrame) * (ps(m, :, iFrame).^2)';
    end
    mu_f = mu_f / sum(ws(:, iFrame));
    sig_f = sig_f / sum(ws(:, iFrame)) - (mu_f.^2);
    mu_fb = mu_fb / sum(wTs(:, iFrame));
    sig_fb = sig_fb / sum(wTs(:, iFrame)) - (mu_fb.^2);
    
    mus_f(:, iFrame) = mu_f;
    sigs_f(:, iFrame) = sqrt(sig_f);
    mus_fb(:, iFrame) = mu_fb;
    sigs_fb(:, iFrame) = sqrt(sig_fb);
end

close(hwb);

close all;

% save('result_Smoother', 'mus_f', 'sigs_f', 'mus_fb', 'sigs_fb');

% Display

vAvg2 = VideoWriter('result_Avg2.avi', 'Uncompressed AVI');
vAvg2.FrameRate = 30;
open(vAvg2);

vTjd2 = VideoWriter('result_Tjd2.avi', 'Uncompressed AVI');
vTjd2.FrameRate = 30;
open(vTjd2);

set(gcf, 'position', [0 0 3*IMWIDTH 2*IMHEIGHT]);
for iFrame=1:nFrame
    subplot(2, 3, 1);
    imshow(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), []);
    subplot(2, 3, 2);
    imshow(impFrame(:, :, iFrame), []);
%     subplot(2, 3, 3);
%     imshow(impBFrame(:, :, iFrame), []);
    subplot(2, 3, 4);
    imshow(mnpFrame(:, :, :, iFrame));
    frm = getframe;
    writeVideo(vAvg2, frm);
    subplot(2, 3, 5);
    imshow(tjd1Frame(:, :, :, iFrame));
    subplot(2, 3, 6);
    imshow(tjd2Frame(:, :, :, iFrame));
    frm = getframe;
    writeVideo(vTjd2, frm);
    
    pause(0.002);
end

close(vAvg2);
close(vTjd2);