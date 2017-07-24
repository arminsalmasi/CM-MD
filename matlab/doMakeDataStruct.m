%% function doMakeDataStruct
% created by Armin Salamsi
% start date 2017-06-09
% 
    % timeSpaceVec.tstep  = timestep #         [tstp         : integer].
    % timeSpaceVec.aIdx   = Index of atoms     [idx          : integer].
    % timeSpaceVec.aTyps  = Atomic number      [atN ...      : integer].
    % timeSpaceVec.aPos   = Pos of Atoms       [x y z ...    : real].
    % timeSpaceVec.aVel   = Velocities.        [vx vy vz ... : real]. 
    % timeSpaceVec.aAcc   = Accelerations.     [ax ay az ... : real].
    % timeSpaceVec.aforce = Interatomic Force. [fx fy fz ... : real].
    % timeSpaceVec.aMass  = Atomic Masses.     [atM ...      : real]. 
    % timeSpaceVec.aRadi  =Atomic radiuses.    [atR ...      : real]. 
    % timeSpaceVec.boxV   = Volume of box.     [boxVol       : real]
    % timeSpaceVec.boxL   = Length of box.     [boxLenght    : real]
    % timeSpaceVec.NP     = # of atoms in box. [Natm         : integer]

function [ a ] = doMakeDataStruct(N,ND, aNums, aXs, T)
    
  a=struct( ...
          'tstep',  [], ...                 %timestep number.  
          'aIdx'  , zeros ( N , 1 ), ...    %Idx of atoms.
          'aTyps' , zeros ( N , 1  ), ...   %Atomic number.
          'aPos'  , zeros ( N , ND ), ...   %Pos of Atoms.
          'aVel'  , zeros ( N , ND ), ...   %Velocities.
          'aAcc'  , zeros ( N , ND ), ...   %Accelerations.
          'aforce', zeros ( N , ND ), ...   %Force.
          'aMass' , zeros ( N , 1  ), ...   %Atomic Masses.
          'aRadi' , zeros ( N , 1  ), ...   %Atomic radiuses.
          'boxV'  , [], ...                 %Volume of box.
          'boxL'  , [], ...                 %Length of box.
          'N'     , N ,...                  %# of atomes
          'ND' , ND, ...                    %# of dimensions.   
          'aNums', aNums, ...               %atomic number of atoms vector
          'aXs', aXs, ...                   % Concentration vector
          'T', T );                         % Prescribed temperature 

end                