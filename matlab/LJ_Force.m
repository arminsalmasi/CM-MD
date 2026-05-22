    function [forces] = LJ_Force(Pos, L, eps, sig)
    
        if nargin < 3
            eps = 1;
        end
        if nargin < 4
            sig = 1;
        end

        isTransposed = false;
        if size(Pos, 1) < size(Pos, 2)
            isTransposed = true;
            Pos = Pos';
        end

        % Initialize all forces to 0
        forces = zeros(size(Pos));
        
        % Get the number of particles
        nPart = size(Pos,1);
        
        sig2 = sig^2;

        % Loop over all particle pairs
        for i = 1:nPart-1
            for j = (i+1):nPart
                
                % Calculate particle-particle distance
                dr = Pos(i,:) - Pos(j,:);
                % Fix according to periodic boundary conditions
                dr = distPBC3D(dr,L);
                % Get the distance squared
                dr2 = dot(dr,dr);
    
                invDr2 = 1.0/dr2; % 1/r^2

                % Performance & Dimensional error fix:
                % F(r) = 24 * eps / r^2 * ( 2*(sig^2/r^2)^6 - (sig^2/r^2)^3 )
                sr2 = sig2 * invDr2;
                sr6 = sr2^3;
                forceFact = 24 * eps * invDr2 * sr6 * (2 * sr6 - 1);
                
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(i,:) = forces(i,:) + dr*forceFact;
                forces(j,:) = forces(j,:) - dr*forceFact;
               
            end
        end
        
        if isTransposed
            forces = forces';
        end
    
    end
