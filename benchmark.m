Pos = rand(500, 3) * 10;
L = 10;
eps = 1;
sig = 1;

tic;
for k = 1:50
    forces = LJ_Force(Pos, L, eps, sig);
end
toc;
