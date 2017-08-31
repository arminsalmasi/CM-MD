% Calculate leonard Jons Potentials
function  [Fcomp] = doPotCalcLJ(a, eps, sig)
  
for i = 1 : a.N-1
  for j= i+1 : a.N
    %% update forces
    if i ~= j
      [rij x y z]= doCalcRij(a, i ,j);
      %U = 4* eps * ( (sig/rij)^12 - (sig/rij)^6 );
      F = ( 24 * eps / sig ) * (2 * (sig/rij)^13 - (sig/rij)^7 ); 
      Fcomp(i,:) =  F * (a.aPos(i,:) - a.aPos(j,:)); Fcomp% [F*(x / rij), F*(y / rij), F*(z / rij)];
      %timeSpaceVec.aforce(i,:) =  Fcomp;
    end
  end
  t=1+1;
end
  
end