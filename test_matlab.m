Pos = rand(100, 3) * 10;
L = 10;
eps = 1;
sig = 1;
disp(sum(sum(LJ_Force(Pos, L, eps, sig))))
