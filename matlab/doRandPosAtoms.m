%% functin dorandPosAtoms - creats a vector which contains randomly
  % distributed atomes based on types of atoms (atomic number) 
% created by Armin Salamsi
% start date 2017-06-12

%% Input
  % timeSpaceVec = first element of datastructure array
  % ND [global] = number of dimensions

%% Output
  % timeSpaceVec = first element of datastructure array

%%
function [timeSpaceVec] = doRandPosAtoms(timeSpaceVec);
  
  %% Randomise stream by system whall clock.
  
  RandStream.setGlobalStream(RandStream('mt19937ar','seed',...
  sum(clock*100))); 
 
  r = zeros ( timeSpaceVec.N, timeSpaceVec.ND );
  a = 0;
  
  %% Loop over all atoms
  for j = 1 : timeSpaceVec.N
    for i = 1 : timeSpaceVec.ND
        r(j,i) = a + (timeSpaceVec.boxL - a) * rand(1);
    end
  end
  %% Report
  timeSpaceVec.aPos = r;
  
end

