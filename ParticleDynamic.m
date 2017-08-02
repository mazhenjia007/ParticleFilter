function [ p_new ] = ParticleDynamic( p_old, M, Height, Width, sigmax, sigmav )
% Execute dynamic process of particles
    p_new = zeros(M, 4);
    
    p_new(:, 1) = p_old(:, 1) + p_old(:, 3) + normrnd(0, sigmax, M, 1);
    p_new(:, 2) = p_old(:, 2) + p_old(:, 4) + normrnd(0, sigmax, M, 1);
    
    p_new(:, 3) = p_old(:, 3) + normrnd(0, sigmav, M, 1);
    p_new(:, 4) = p_old(:, 4) + normrnd(0, sigmav, M, 1);
    
    for m=1:M
        if (p_new(m, 1) <= 0)
            p_new(m, 1) = 0 + abs(normrnd(0, sigmax));
            
            if (p_new(m, 3) < 0)
                p_new(m, 3) = -p_new(m, 3);
            end
        end
        
        if (p_new(m, 1) >= Height)
            p_new(m, 1) = Height - abs(normrnd(0, sigmax));
            
            if (p_new(m, 3) > 0)
                p_new(m, 3) = -p_new(m, 3);
            end
        end

        if (p_new(m, 2) <= 0)
            p_new(m, 2) = 0 + abs(normrnd(0, sigmax));
            
            if (p_new(m, 4) < 0)
                p_new(m, 4) = -p_new(m, 4);
            end
        end
        
        if (p_new(m, 2) >= Width)
            p_new(m, 2) = Width - abs(normrnd(0, sigmax));
            
            if (p_new(m, 4) > 0)
                p_new(m, 4) = -p_new(m, 4);
            end
        end
    end
end

