clear variables
close all
clc
coords = [1 2 3 4 5 ; 1 2 3 4 5 ; 1 2 3 4 5]
L = 1 
forces = LJ_Force(coords,L)

function forces = LJ_Force(coords,L)
      
        % Initialize all forces to 0
        forces = zeros(size(coords))
        
        % Get the number of particles
        nPart = size(coords,2)
        k = 1
        % Loop over all particle pairs
        for partA = 1:nPart-1
            for partB = (partA+1):nPart
                partA
                partB
                % Calculate particle-particle distance
                dr = coords(:,partA) - coords(:,partB)
                % Fix according to periodic boundary conditions
                %dr = distPBC3D(dr,L);
                % Get the distance squared
                dr2 = dot(dr,dr)
    
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
                    
                invDr2 = 1.0/dr2 % 1/r^2
                forceFact = invDr2^4 * (invDr2^3 - 0.5)
                    
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(:,partA) = forces(:,partA) + dr*forceFact
                forces(:,partB) = forces(:,partB) - dr*forceFact
               k = k+1
            end
        end
        
        % Multiply all forces by 48
        forces = forces*48;
    
    end