%% functin doFixCoM - fixes center of mass to prevent drift
  %Set initial momentum to zero.
% created by Armin Salamsi
% start date 2017-06-12

% Calcualte velocity of the center of mass VCM = sum(mi*Vi)/sum(mi) .
% Reduce the velocity to fix the center of mass    

%% Input
  % timeSpaceVec = first element of datastructure array
  % T [global] = temperature [K]
  % ND [global] = Number of Dimentions
  % KB [global] = Boltzmann constant [m^2.kg.s^-2.K^-1].  
%% Output
  % timeSpaceVec = first element of datastructure array

%%
function [timeSpaceVec] = doFixCoM(timeSpaceVec)
  %% Calcualte velocity of the center of mass VCM = sum(mi*Vi)/sum(mi) .
  %% loop over all atoms
    for i = 1 :  timeSpaceVec.N 
      vCoMtemp(i,:) =  timeSpaceVec.aVel(i,:) .* timeSpaceVec.aMass(i) ; 
    end
    
  %% calculate velocity of center of mass
    vCoM = sum(vCoMtemp) / sum(timeSpaceVec.aMass);
    
  %% Fix any center-of-mass drift V= V - VCoM  
  %% Loop over all atoms and fix the center of mass
    for i= 1 : timeSpaceVec.N                   %    
      timeSpaceVec.aVel(i,:) = timeSpaceVec.aVel(i,:) - vCoM;    
    end
    
end

