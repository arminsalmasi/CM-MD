%% functin doScaleVel - scale velocities by 
  %the prescribed temperature factor = T/T0
% created by Armin Salamsi
% start date 2017-06-12

% Set initial kinetic energy to nd*kBT/2.
% sum(dot(v,v))= v(i,1)^2+v(i,2)^2+v(i,3)^2.  = V mean squered

%% Input
  % timeSpaceVec = first element of datastructure array
  % T [global] = temperature [K]
  % ND [global] = Number of Dimentions
  % KB [global] = Boltzmann constant [m^2.kg.s^-2.K^-1].  

%% Output
  % timeSpaceVec = first element of datastructure array

%%
function [timeSpaceVec] = doScaleVel(timeSpaceVec)
  
  
  n=1; %counter
  %% calculate the temperature 
    Ttemp = doTCalc(timeSpaceVec);
  
  %% correct the velocities to fit the temperature (itterations)
    TCorrectionFactor = ((timeSpaceVec.T -Ttemp) /timeSpaceVec.T) *100; 
    %% itterates to reach the prescribed accuracy or number of maximum 
      % number itteration
    while ( (abs(TCorrectionFactor)>0.1) && (n <1000) )
      vScale = sqrt(timeSpaceVec.T/Ttemp);
      for i = 1:timeSpaceVec.N                          
        timeSpaceVec.aVel(i,:) = timeSpaceVec.aVel(i,:)* vScale;    
      end
      Ttemp = doTCalc(timeSpaceVec);
      TCorrectionFactor = ((timeSpaceVec.T-Ttemp)/timeSpaceVec.T)*100;
      n = n+1;
  end
  
end

