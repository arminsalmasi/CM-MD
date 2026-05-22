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

  %% Randomise stream securely.
  % Security fix: Use java.security.SecureRandom instead of predictable sum(clock*100)
  % to prevent DoS or simulation manipulation via predictable seed values.
  try
      secureRandom = java.security.SecureRandom();
      % Generate 4 secure bytes and cast to uint32 for the seed
      seedBytes = secureRandom.generateSeed(4);
      seed = typecast(seedBytes, 'uint32');
      RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', seed));
  catch
      % Fallback to shuffle if Java is unavailable
      rng('shuffle', 'twister');
  end

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
