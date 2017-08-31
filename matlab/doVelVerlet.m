function [Pos, Vel, Acc]= doVelVerlet(Pos, Vel, Acc, force, Mass, dt);
   
  Pos = Pos + Vel * dt + 0.5 * Acc * dt * dt;

  for i = 1: size(Mass,1)
    Vel(i,:) = Vel(i,:) + 0.5 * ( force(i,:) / Mass(i) + Acc(i,:)) *dt;
  end
  
  for i = 1 : size(Mass,1)
    Acc(i,:) =  force(i,:) / Mass(i);
  end
end