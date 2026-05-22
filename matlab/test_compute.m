function test_compute()
    % Basic test case
    pos1 = [0, 1; 0, 0; 0, 0];
    vel1 = [1, -1; 0, 0; 0, 0];
    mass1 = 2.0;
    [f1, pot1, kin1] = compute(2, 3, pos1, vel1, mass1);

    assert(abs(pot1 - 0.7080734182735712) < 1e-6, 'Basic potential energy failed');
    assert(abs(kin1 - 2.0) < 1e-6, 'Basic kinetic energy failed');
    assert(abs(f1(1,1) - 0.9092974268256817) < 1e-6, 'Basic force X1 failed');
    assert(abs(f1(1,2) - -0.9092974268256817) < 1e-6, 'Basic force X2 failed');
    assert(abs(f1(2,1) - 0.0) < 1e-6, 'Basic force Y1 failed');
    assert(abs(f1(3,1) - 0.0) < 1e-6, 'Basic force Z1 failed');

    % Far distance test case
    pos2 = [0, 2; 0, 0; 0, 0];
    [f2, pot2, kin2] = compute(2, 3, pos2, vel1, mass1);

    assert(abs(pot2 - 1.0) < 1e-6, 'Far potential energy failed');
    assert(abs(f2(1,1) - 0.0) < 1e-6, 'Far force X1 failed');

    disp('All tests passed.');
end
