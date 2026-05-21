% Calculate leonard Jons Potentials
function  [Fcomp] = doPotCalcLJ(a, eps, sig)

for i = 1 : a.N-1
  for j= i+1 : a.N
    %% update forces
    if i ~= j
      [rij x y z]= doCalcRij(a, i ,j);
      %U = 4* eps * ( (sig/rij)^12 - (sig/rij)^6 );

      % Performance optimization: avoid large exponentiation
      % by fixing dimensional error to use squared distance.
      rij2 = rij^2;
      invRij2 = 1.0 / rij2;
      sig2InvRij2 = (sig^2) * invRij2;
      sr6 = sig2InvRij2 * sig2InvRij2 * sig2InvRij2;
      F_r_over_r = (48 * eps) * invRij2 * sr6 * (sr6 - 0.5);
      Fcomp(i,:) =  F_r_over_r * (a.aPos(i,:) - a.aPos(j,:)); Fcomp% [F*(x / rij), F*(y / rij), F*(z / rij)];
      %timeSpaceVec.aforce(i,:) =  Fcomp;
    end
  end
  t=1+1;
end

end