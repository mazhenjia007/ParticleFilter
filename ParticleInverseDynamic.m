function [ p_old ] = ParticleInverseDynamic( p_new, M, Height, Width, sigmax, sigmav, p_f, w_f )
% Execute inverse dynamic process of particles

    A = [   1, 0, 1, 0;
            0, 1, 0, 1;
            0, 0, 1, 0;
            0, 0, 0, 1];
        
    Q = diag([sigmax^2, sigmax^2, sigmav^2, sigmav^2]);
    Mu_f = zeros(4, 1);
    Sigma_f = zeros(4, 4);
    
    for m=1:M
        Mu_f = Mu_f + w_f(m)*p_f(m, :)';
        Sigma_f = Sigma_f + w_f(m)*p_f(m, :)'*p_f(m, :);
    end
    Mu_f = Mu_f / sum(w_f);
    Sigma_f = Sigma_f / sum(w_f) - Mu_f*Mu_f';
    
    p_old = zeros(size(p_new));
    iSigma = inv(Sigma_f) + A'/Q*A;
    Sigma = inv(iSigma);
    Sigma = (Sigma + Sigma') /2;
    for m=1:M
        Mu = iSigma \ (A'/Q*p_new(m, :)' + Sigma_f\Mu_f);
        % Mu = Mu + p_new(m, :)';
        p_old(m, :) = mvnrnd(Mu, Sigma);
    end
    
    for m=1:M
        if (p_old(m, 1) <= 0)
            p_old(m, 1) = 0 + abs(normrnd(0, sigmax));
            
            if (p_old(m, 3) < 0)
                p_old(m, 3) = -p_old(m, 3);
            end
        end
        
        if (p_old(m, 1) >= Height)
            p_old(m, 1) = Height - abs(normrnd(0, sigmax));
            
            if (p_old(m, 3) > 0)
                p_old(m, 3) = -p_old(m, 3);
            end
        end

        if (p_old(m, 2) <= 0)
            p_old(m, 2) = 0 + abs(normrnd(0, sigmax));
            
            if (p_old(m, 4) < 0)
                p_old(m, 4) = -p_old(m, 4);
            end
        end
        
        if (p_old(m, 2) >= Width)
            p_old(m, 2) = Width - abs(normrnd(0, sigmax));
            
            if (p_old(m, 4) > 0)
                p_old(m, 4) = -p_old(m, 4);
            end
        end
    end
end

