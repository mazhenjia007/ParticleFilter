load('img400b.mat');

IMHEIGHT = 300;
IMWIDTH = 400;

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

ds = 1;
iss = ds:ds:Height;
jss = ds:ds:Width;
niss = Height / ds;
njss = Width / ds;
llhFrame = zeros(niss, njss, nFrame);

% #particles
M = 5000;

% [x1, x2, v1, v2] : x1 v1 in Height direction, x2 v2 in Width direction

sigmax = 0.05;
sigmav = 0.1;

p = zeros(M, 4);
w = zeros(M, 1);

p(:, 1) = rand(M, 1) * Height;
p(:, 2) = rand(M, 1) * Width;

% particles = SetInitParticles(M, Height, Width);

% figId = figure(1);

impFrame = zeros(IMHEIGHT, IMWIDTH, nFrame);

mnpFrame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);

mleFrame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);

tjdFrame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
lastFrame = zeros(IMHEIGHT, IMWIDTH, 3);

hwb = waitbar(0, 'Initializing...');

tag = 0;

lastx = -1;
lasty = -1;

for iFrame = 1:nFrame

    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('%d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));

    % Dynamic   
    p = ParticleDynamic(p, M, Height, Width, sigmax, sigmav);
        
    % Weighting ( Likelihood )
    for m=1:M
        w(m) = CalLikelihood_Subpixel(vs(:, :, iFrame), ...
            Rad, p(m, 1), p(m, 2));
        w(m) = exp(w(m));
    end
    
    % Resampling
    plabel = randsample(1:M, M, true, w);
    p = p(plabel, :);
    
    % Whole Likelihood
    llhi =zeros(niss, njss);
    for i=1:niss
        for j=1:njss
            
            lln = CalLikelihood_Subpixel(vs(:, :, iFrame), ...
                Rad, iss(i), jss(j));
            llhi(i, j) = exp(lln);
        end
    end
    llhFrame(:, :, iFrame) = llhi / sum(sum(llhi));

    % Display
%     subplot(1, 2, 1);
%     imshow(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), []);

    imp = zeros(IMHEIGHT, IMWIDTH);
    for m=1:M
        imp(ceil(p(m, 1)/Height*IMHEIGHT), ceil(p(m, 2)/Width*IMWIDTH)) = 1;
    end
    impFrame(:, :, iFrame) = imp;
    
    mnp = ind2rgb(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), [0, 0, 0; 1, 1, 1]);
    mx = sum(p(:, 1).*w) / sum(w) / Height * IMHEIGHT;
    my = sum(p(:, 2).*w) / sum(w) / Width * IMWIDTH;
    mnp = DrawBox(mnp, IMHEIGHT, IMWIDTH, mx, my, Rad/Height*IMHEIGHT, Rad/Width*IMWIDTH);
    mnpFrame(:, :, :, iFrame) = mnp;
    
    if (lastx ~= -1)
        dx = mx - lastx;
        dy = my - lasty;
        dt = min(floor(abs(dx/Height*IMHEIGHT)), floor(abs(dy/Width*IMWIDTH)));
        dx = dx / dt;
        dy = dy / dt;
        for t=1:dt
            lastFrame = DrawBox(lastFrame, IMHEIGHT, IMWIDTH, lastx+dx*t, lasty+dy*t, 0, 0);
        end
    end
    lastFrame = DrawBox(lastFrame, IMHEIGHT, IMWIDTH, mx, my, 0, 0);
    tjdFrame(:, :, :, iFrame) = lastFrame;
    lastx = mx; lasty = my;
    
    mle = ind2rgb(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), [0, 0, 0; 1, 1, 1]);
    [~, indp] = sort(w);
    
    mx = p(indp(1), 1) / Height * IMHEIGHT;
    my = p(indp(1), 2) / Width * IMWIDTH;
    mle = DrawBox(mle, IMHEIGHT, IMWIDTH, mx, my, Rad/Height*IMHEIGHT, Rad/Width*IMWIDTH);
    mx = p(indp(2), 1) / Height * IMHEIGHT;
    my = p(indp(2), 2) / Width * IMWIDTH;
    mle = DrawBox(mle, IMHEIGHT, IMWIDTH, mx, my, Rad/Height*IMHEIGHT, Rad/Width*IMWIDTH);
    mx = p(indp(3), 1) / Height * IMHEIGHT;
    my = p(indp(3), 2) / Width * IMWIDTH;
    mle = DrawBox(mle, IMHEIGHT, IMWIDTH, mx, my, Rad/Height*IMHEIGHT, Rad/Width*IMWIDTH);
    mleFrame(:, :, :, iFrame) = mle;
%     subplot(1, 2, 2);
%     imshow(imp, []);
%     pause(0.001);
end

close(hwb);

save('result', 'impFrame', 'llhFrame', 'tjdFrame');

close all;

% Display

vAvg1 = VideoWriter('result_Avg1.avi', 'Uncompressed AVI');
vAvg1.FrameRate = 30;
open(vAvg1);

vMLE = VideoWriter('result_MLE.avi', 'Uncompressed AVI');
vMLE.FrameRate = 30;
open(vMLE);

vTjd1 = VideoWriter('result_Tjd1.avi', 'Uncompressed AVI');
vTjd1.FrameRate = 30;
open(vTjd1);

set(gcf, 'position', [0 0 2*IMWIDTH 3*IMHEIGHT]);
for iFrame=1:nFrame
    subplot(3, 2, 1);
    imshow(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), []);
    subplot(3, 2, 2);
    imshow(imresize(llhFrame(:, :, iFrame), [IMHEIGHT, IMWIDTH]), []);
    subplot(3, 2, 3);
    imshow(impFrame(:, :, iFrame), []);
    subplot(3, 2, 4);
    imshow(mnpFrame(:, :, :, iFrame));
    frm = getframe;
    writeVideo(vAvg1, frm);
    subplot(3, 2, 5);
    imshow(mleFrame(:, :, :, iFrame));
    frm = getframe;
    writeVideo(vMLE, frm);
    subplot(3, 2, 6);
    imshow(tjdFrame(:, :, :, iFrame));
    frm = getframe;
    writeVideo(vTjd1, frm);
    
    pause(0.002);
end

close(vAvg1);
close(vMLE);
close(vTjd1);