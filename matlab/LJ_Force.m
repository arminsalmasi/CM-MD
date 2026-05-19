    function [forces] = LJ_Force(Pos,L, eps, sig)
    
        % Initialize all forces to 0
        forces = zeros(size(Pos));
        
        % Get the number of particles
        nPart = size(Pos,1);
        
        % Precompute sigma squared
        sig2 = sig^2;

        % Loop over all particle pairs
        for i = 1:nPart-1
            for j = (i+1):nPart
                
                % Calculate particle-particle distance
                dr = Pos(i,:) - Pos(j,:);
                % Fix according to periodic boundary conditions
                dr = distPBC3D(dr,L);

                % Get the distance squared directly (avoiding sqrt)
                dr2 = dr(1)^2 + dr(2)^2 + dr(3)^2;
    
                % Lennard-Jones potential:
                % U(r) = 4*epsilon* [(sigma/r)^12 - (sigma/r)^6]
                % Fx(r) = 48 * epsilon * x * (1/r^2) * [ (sigma/r)^12 - 0.5 * (sigma/r)^6 ]
                %
                % For efficiency, we avoid explicit exponentiation and sqrt,
                % and multiply by 48 * epsilon only after summing up all the forces.
                    
                invDr2 = 1.0 / dr2; % 1/r^2
                sr2 = sig2 * invDr2; % (sigma/r)^2
                sr6 = sr2 * sr2 * sr2; % (sigma/r)^6
                forceFact = invDr2 * sr6 * (sr6 - 0.5);
                
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(i,:) = forces(i,:) + dr*forceFact;
                forces(j,:) = forces(j,:) - dr*forceFact;
               
            end
        end
        
        % Multiply all forces by 48 * epsilon
        forces = forces * 48 * eps;
    
    end
