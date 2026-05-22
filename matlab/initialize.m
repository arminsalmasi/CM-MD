function [ pos, vel, acc ] = initialize ( np, nd )

%*****************************************************************************80
%
%% INITIALIZE initializes the positions, velocities, and accelerations.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 December 2014
%
%  Author:
%
%    John Burkardt.
%
%  Parameters:
%
%    Input, integer NP, the number of particles.
%
%    Input, integer ND, the number of spatial dimensions.
%
%    Output, real POS(ND,NP), the positions.
%
%    Output, real VEL(ND,NP), the velocities.
%
%    Output, real ACC(ND,NP), the accelerations.
%
  seed = 123456789;
%
%  Positions.
%
  [ pos, seed ] = r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed );
%
%  Velocities.
%
  vel = zeros ( nd, np );
%
%  Accelerations.
%
  acc = zeros ( nd, np );

  return
end
