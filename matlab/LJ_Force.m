    function [forces] = LJ_Force(Pos,L, eps, sig)
    
        % Initialize all forces to 0
        forces = zeros(size(Pos));
        
        % Get the number of particles
        nPart = size(Pos,1);
        
        % Precalculate constants for performance
        hL = L/2.0;
        sig2 = sig*sig;

        % Loop over all particle pairs
        for i = 1:nPart-1
            for j = (i+1):nPart
                
                % Calculate particle-particle distance
                dr = Pos(i,:) - Pos(j,:);

                % Fix according to periodic boundary conditions (inlined)
                if dr(1) > hL, dr(1) = dr(1) - L; elseif dr(1) < -hL, dr(1) = dr(1) + L; end
                if dr(2) > hL, dr(2) = dr(2) - L; elseif dr(2) < -hL, dr(2) = dr(2) + L; end
                if dr(3) > hL, dr(3) = dr(3) - L; elseif dr(3) < -hL, dr(3) = dr(3) + L; end

                % Get the distance squared
                dr2 = dr(1)^2 + dr(2)^2 + dr(3)^2;
    
                % Lennard-Jones potential:
                % U(r) = 4*epsilon* [(sigma/r)^12 - (sigma/r)^6]
                %
                % Here, we set sigma = 1, epsilon = 1 (reduced distance and
                % energy units). Therefore:
                %
                % U(r) = 4 * [(1/r)^12 - (1/r)^6]
                % 
                % Fx(r) = 4 * (x/r) * [12*(1/r)^13 - 6*(1/r)^7]
                %
                %      = 48 * x * (1/r)^8 * [(1/r)^6 - 0.5]
                %
                % And same goes for the force in the y&z directions.
                %
                % For efficiency, we will multiply by 48 only after summing
                % up all the forces.
                    
                %invDr2 = 1.0/dr2; % 1/r^2
                %forceFact = invDr2^4 * (invDr2^3 - 0.5);

                % Performance optimization: avoid large exponents
                invDr2 = sig2 / dr2;
                invDr6 = invDr2 * invDr2 * invDr2;
                forceFact = 2 * sqrt(invDr2) * invDr6 * (invDr6 - 1.0);
                
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(i,:) = forces(i,:) + dr*forceFact;
                forces(j,:) = forces(j,:) - dr*forceFact;
               
            end
        end
        
        % Multiply all forces by 48
        forces = forces*48 * eps /sig;
    
    end