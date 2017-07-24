    function [forces] = LJ_Force(Pos,L, eps, sig)
    
        % Initialize all forces to 0
        forces = zeros(size(Pos));
        
        % Get the number of particles
        nPart = size(Pos,1);
        
        % Loop over all particle pairs
        for i = 1:nPart-1
            for j = (i+1):nPart
                
                % Calculate particle-particle distance
                dr = Pos(i,:) - Pos(j,:);
                % Fix according to periodic boundary conditions
                dr = distPBC3D(dr,L);
                % Get the distance squared
                %dr2 = dot(dr,dr);
    
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
                drsize = sqrt(dr(1)^2+dr(2)^2+dr(3)^2);
                    forceFact = 2 * ( (sig/drsize)^13 - (sig/drsize)^7 );    
                
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(i,:) = forces(i,:) + dr*forceFact;
                forces(j,:) = forces(j,:) - dr*forceFact;
               
            end
        end
        
        % Multiply all forces by 48
        forces = forces*48 * eps /sig;
    
    end