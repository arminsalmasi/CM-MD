%% functin doRandVel - cprescribe a random velocity vector based on the
  %Temperature 
% created by Armin Salamsi
% start date 2017-06-12

%% Input
  % timeSpaceVec = first element of datastructure array
  % T [global] = temperature [K]
  % ND [global] = Number of Dimentions
  % KB [global] = Boltzmann constant [m^2.kg.s^-2.K^-1].  

%% Output
  % timeSpaceVec = first element of datastructure array


%%  
function [timeSpaceVec] = doRandVel(timeSpaceVec)
  global  KB 
  
  %% loop over all atoms  
  for  i = 1 : timeSpaceVec.N
    % vSTD = standard deviation of velocities of all atoms
    vSTD = sqrt( timeSpaceVec.T * KB  / timeSpaceVec.aMass( i ) );
    timeSpaceVec.aVel(i,:) = 0 + randn( 1 , timeSpaceVec.ND) * vSTD; 
  end
  
end



