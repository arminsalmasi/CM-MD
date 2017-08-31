%% functin doTCalc - Calculates temperatures

% created by Armin Salamsi
% start date 2017-06-12

% T = sum( Mass(i) *V(i)) /(KB *Nf) ) 

%% Input
  % a = timeSpaceVec = first element of datastructure array

%% Output
  % calculated temperature after correction itterations

%%
function [T] = doTCalc(a)
  global KB 
  
  % Nf = number of degrees of freedom
  Nf = 3 *a.N - 3;
  
  % starting value of T summation
  T = 0;
  
  %% loop over all atoms
  for i = 1 : a.N
    % magnitude of the velocity vectores ( all atomes on by one)   
    ViSq =  sum( dot( a.aVel(i ,:) ,a.aVel(i ,:) ) ); 
    T = T +(a.aMass(i) *ViSq) /(KB *Nf);
    ViSq = 0;
  end
  
end