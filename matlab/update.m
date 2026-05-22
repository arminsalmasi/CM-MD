function [ pos, vel, acc ] = update ( np, nd, pos, vel, f, acc, mass, dt )

%*****************************************************************************80
%
%% UPDATE updates positions, velocities and accelerations.
%
%  Discussion:
%
%    The time integration is fully parallel.
%
%    A velocity Verlet algorithm is used for the updating.
%
%    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
%    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
%    a(t+dt) = f(t) / m
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    15 July 2008
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NP, the number of particles.
%
%    Input, integer ND, the number of spatial dimensions.
%
%    Input, real POS(ND,NP), the positions.
%
%    Input, real VEL(ND,NP), the velocities.
%
%    Input, real F(ND,NP), the forces.
%
%    Input, real ACC(ND,NP), the accelerations.
%
%    Input, real MASS, the mass of each particle.
%
%    Input, real DT, the time step.
%
%    Output, real POS(ND,NP), the updated positions.
%
%    Output, real VEL(ND,NP), the updated velocities.
%
%    Output, real ACC(ND,NP), the updated accelerations.
%
  rmass = 1.0 / mass;

  pos(1:nd,1:np) = pos(1:nd,1:np) + vel(1:nd,1:np) * dt ...
    + 0.5 * acc(1:nd,1:np) * dt * dt;

  vel(1:nd,1:np) = vel(1:nd,1:np) ...
    + 0.5 * dt * ( f(1:nd,1:np) * rmass + acc(1:nd,1:np) );

  acc(1:nd,1:np) = f(1:nd,1:np) * rmass;

  return
end