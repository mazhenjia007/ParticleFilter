load('img400b.mat');

IMHEIGHT = 300;
IMWIDTH = 400;

eps = 1e-6;

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
M = 1000;

% [x1, x2, v1, v2] : x1 v1 in Height direction, x2 v2 in Width direction

sigmax = 0.05;
sigmav = 0.1;

p = zeros(M, 4);
w = zeros(M, 1);

p(:, 1) = rand(M, 1) * Height;
p(:, 2) = rand(M, 1) * Width;

p0 = p;
w0 = ones(M, 1) / M;

ps = zeros(M, 4, nFrame);
ws = zeros(M, nFrame);
p_backup = zeros(M, 4);
w_backup = zeros(M, 1);

% particles = SetInitParticles(M, Height, Width);

% figId = figure(1);

hwb = waitbar(0, 'Initializing...');

tag = 0;

for iFrame = 1:nFrame
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Forward: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));

    if (iFrame == nFrame)
        p_backup = p;
    end
        
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
    
    if (iFrame == nFrame-1)
        w_backup = w(plabel, :);
    end
end

Mus = zeros(4, nFrame);
Sigs = zeros(4, 4, nFrame);
for iFrame=1:nFrame
    Mu = zeros(4, 1);
    Sig = zeros(4, 4);
    for n=1:M
        Mu = Mu + ws(n, iFrame)*ps(n, :, iFrame)';
        Sig = Sig + ws(n, iFrame)*ps(n, :, iFrame)'*ps(n, :, iFrame);
    end
    Mu = Mu / sum(ws(:, iFrame));
    Sig = Sig / sum(ws(:, iFrame)) - Mu*Mu';
    
    Sig = (Sig + Sig') / 2;
    
    Mus(:, iFrame) = Mu;
    Sigs(:, :, iFrame) = Sig;
end

% wT = zeros(M, 1);
wTs = zeros(M, nFrame);

wT = w;

p = ps(:, :, nFrame);
for m=1:M
    pox = p_backup(m, 1:2) + p_backup(m, 3:4);
    pov = p_backup(m, 3:4);
    pnx = p(m, 1:2);
    pnv = p(m, 3:4);
    
    dx = (pnx - pox) / sigmax;
    dv = (pnv - pov) / sigmav;
    
    dlt = (- dx(:, 1).^2 - dx(:, 2).^2 - dv(:, 1).^2 - dv(:, 2).^2) / 2;
    
    if (exp(dlt) < eps)
        wT(m) = eps;
    else
        wT(m) = ws(m, nFrame) / w_backup(m) / exp(dlt);
    end
end

wT = wT / sum(wT);

wTs(:, nFrame) = wT;

pBs = zeros(M, 4, nFrame);
pBs(:, :, nFrame) = ps(:, :, nFrame);

tag = 0;
for iFrame=nFrame-1:-1:1
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Backward: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));
        
    % display(p);
    
    p_backup = p;
    
    % Inverse Dynamic   
    p = ParticleInverseDynamic(p, M, Height, Width, sigmax, sigmav, ps(:, :, iFrame), ws(:, iFrame));
        
    % Weighting ( Likelihood )
    for m=1:M
        wT(m) = CalLikelihood_Subpixel(vs(:, :, iFrame), ...
            Rad, p(m, 1), p(m, 2));
        wT(m) = exp(wT(m));
        
%         pdf1 = mvnpdf(p(m, :)', Mus(:, iFrame), Sigs(:, :, iFrame));
%         pdf2 = mvnpdf(p_backup(m, :)', Mus(:, iFrame+1), Sigs(:, :, iFrame+1));
%         wT(m) = wT(m) * pdf1 / pdf2; 
    end
    
    wT = wT / sum(wT);
    
    pBs(:, :, iFrame) = p;
    wTs(:, iFrame) = wT;
    
    % Resampling
    plabel = randsample(1:M, M, true, wT);
    p = p(plabel, :);
end

nu = zeros(M, 1);
nus = zeros(M, nFrame);

tag = 0;
for iFrame=1:nFrame
    tag = tag + 1;
    perc = tag / (nFrame);
    waitbar(perc, hwb, sprintf('Forward-Backward: %d/%d, %.2f%% ...', tag,...
            nFrame, perc*100));
        
    Mu = Mus(:, iFrame);
    Sig = Sigs(:, :, iFrame);
    
    for m=1:M
        pmx = pBs(m, 1:2, iFrame);
        pmv = pBs(m, 3:4, iFrame);
        
        if (iFrame == 1)
            pnx = p0(:, 1:2) + p0(:, 3:4);
            pnv = p0(:, 3:4);
            wn = w0;
        else
            pnx = ps(:, 1:2, iFrame-1) + ps(:, 3:4, iFrame-1);
            pnv = ps(:, 3:4, iFrame-1);
            wn = ws(:, iFrame-1);
        end
        
        dx = bsxfun(@minus, pmx, pnx) / sigmax;
        dv = bsxfun(@minus, pmv, pnv) / sigmav;
        
        dlt = (- dx(:, 1).^2 - dx(:, 2).^2 - dv(:, 1).^2 - dv(:, 2).^2) / 2;
        % dlt = (- dx(:, 1).^2 - dx(:, 2).^2) / 2;
        pdf1 = bsxfun(@times, wn, exp(dlt));
        pdf1 = sum(sum(pdf1));
        
        pdf2 = mvnpdf(pBs(m, :, iFrame)', Mu, Sig);
        % pdf2 = mvnpdf(pBs(m, 1:2, iFrame)', Mu(1:2), Sig(1:2, 1:2));
        
        if (pdf2 < eps)
            nu(m) = 0;
        else
            if (pdf1 < eps)
                nu(m) = 1e-8;
            else
                nu(m) = wTs(m, iFrame) * pdf1 / pdf2;
            end
        end
    end
    
    nu = nu / sum(sum(nu));
    nus(:, iFrame) = nu;
end

impFrame = zeros(IMHEIGHT, IMWIDTH, nFrame);
impBFrame = zeros(IMHEIGHT, IMWIDTH, nFrame);

lastx1 = -1;
lasty1 = -1;
lastx2 = -1;
lasty2 = -1;
lastx3 = -1;
lasty3 = -1;
tjd1Frame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
tjd2Frame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
tjd3Frame = zeros(IMHEIGHT, IMWIDTH, 3, nFrame);
lastFrame1 = zeros(IMHEIGHT, IMWIDTH, 3);
lastFrame2 = zeros(IMHEIGHT, IMWIDTH, 3);
lastFrame3 = zeros(IMHEIGHT, IMWIDTH, 3);

Mus_f = zeros(4, nFrame);
Sigs_f = zeros(4, nFrame);
Mus_b = zeros(4, nFrame);
Sigs_b = zeros(4, nFrame);
Mus_fb = zeros(4, nFrame);
Sigs_fb = zeros(4, nFrame);

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
    
    impB = zeros(IMHEIGHT, IMWIDTH);
    for m=1:M
        impB(ceil(pBs(m, 1, iFrame)/Height*IMHEIGHT), ceil(pBs(m, 2, iFrame)/Width*IMWIDTH)) = 1;
    end
    impBFrame(:, :, iFrame) = impB;
    
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
    
    mx2 = sum(pBs(:, 1, iFrame).*wTs(:, iFrame)) / sum(wTs(:, iFrame)) / Height * IMHEIGHT;
    my2 = sum(pBs(:, 2, iFrame).*wTs(:, iFrame)) / sum(wTs(:, iFrame)) / Width * IMWIDTH;

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

    mx3 = sum(pBs(:, 1, iFrame).*nus(:, iFrame)) / sum(nus(:, iFrame)) / Height * IMHEIGHT;
    my3 = sum(pBs(:, 2, iFrame).*nus(:, iFrame)) / sum(nus(:, iFrame)) / Width * IMWIDTH;

    if (lastx3 ~= -1)
        dx = mx3 - lastx3;
        dy = my3 - lasty3;
        dt = min(floor(abs(dx/Height*IMHEIGHT)), floor(abs(dy/Width*IMWIDTH)));
        dx = dx / dt;
        dy = dy / dt;
        for t=1:dt
            lastFrame3 = DrawBox(lastFrame3, IMHEIGHT, IMWIDTH, lastx3+dx*t, lasty3+dy*t, 0, 0);
        end
    end
    lastFrame3 = DrawBox(lastFrame3, IMHEIGHT, IMWIDTH, mx3, my3, 0, 0);
    tjd3Frame(:, :, :, iFrame) = lastFrame3;
    lastx3 = mx3;
    lasty3 = my3;

    Mu_f = zeros(4, 1);
    Sig_f = zeros(4, 1);
    Mu_b = zeros(4, 1);
    Sig_b = zeros(4, 1);
    Mu_fb = zeros(4, 1);
    Sig_fb = zeros(4, 1);
    for m=1:M
        Mu_f = Mu_f + ws(m, iFrame) * ps(m, :, iFrame)';
        Sig_f = Sig_f + ws(m, iFrame) * (ps(m, :, iFrame).^2)';
        Mu_b = Mu_b + wTs(m, iFrame) * pBs(m, :, iFrame)';
        Sig_b = Sig_b + wTs(m, iFrame) * (pBs(m, :, iFrame).^2)';
%         Mu_fb = Mu_fb + nus(m, iFrame) * pBs(m, :, iFrame)';
%         Sig_fb = Sig_fb + nus(m, iFrame) * (pBs(m, :, iFrame).^2)';
    end
    Mu_f = Mu_f / sum(ws(:, iFrame));
    Sig_f = Sig_f / sum(ws(:, iFrame)) - (Mu_f.^2);
    Mu_b = Mu_b / sum(wTs(:, iFrame));
    Sig_b = Sig_b / sum(wTs(:, iFrame)) - (Mu_b.^2);
%     Mu_fb = Mu_fb / sum(nus(:, iFrame));
%     Sig_fb = Sig_fb / sum(nus(:, iFrame)) - (Mu_fb.^2);
    
    Mus_f(:, iFrame) = Mu_f;
    Sigs_f(:, iFrame) = Sig_f;
    Mus_b(:, iFrame) = Mu_b;
    Sigs_b(:, iFrame) = Sig_b;    
    Mus_fb(:, iFrame) = Mu_fb;
    Sigs_fb(:, iFrame) = Sig_fb;    
end

close(hwb);

save('result_Smoother2', 'Mus_f', 'Sigs_f', 'Mus_b', 'Sigs_b', 'Mus_fb', 'Sigs_fb');

close all;

% Display

set(gcf, 'position', [0 0 3*IMWIDTH 2*IMHEIGHT]);
for iFrame=1:nFrame
    subplot(2, 3, 1);
    imshow(imresize(vs(:, :, iFrame), [IMHEIGHT, IMWIDTH]), []);
    subplot(2, 3, 2);
    imshow(impFrame(:, :, iFrame), []);
    subplot(2, 3, 3);
    imshow(impBFrame(:, :, iFrame), []);
    subplot(2, 3, 4);
    imshow(tjd3Frame(:, :, :, iFrame));
    subplot(2, 3, 5);
    imshow(tjd1Frame(:, :, :, iFrame));
    subplot(2, 3, 6);
    imshow(tjd2Frame(:, :, :, iFrame));
    pause(0.002);
end
