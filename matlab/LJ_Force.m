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
                dr2 = dr(1)^2 + dr(2)^2 + dr(3)^2;

                % Performance optimization: avoid sqrt() and large exponentiation
                % by fixing dimensional error to use squared distance.
                invDr2 = 1.0/dr2; % 1/r^2
                sig2InvDr2 = (sig^2)*invDr2; % (sigma/r)^2
                sr6 = sig2InvDr2 * sig2InvDr2 * sig2InvDr2; % (sigma/r)^6
                forceFact = invDr2 * sr6 * (sr6 - 0.5);
                
                % According to Newton's third law, we get action and
                % reaction for the two particles.
                forces(i,:) = forces(i,:) + dr*forceFact;
                forces(j,:) = forces(j,:) - dr*forceFact;
               
            end
        end
        
        % Multiply all forces by 48
        forces = forces * 48 * eps;
    
    end