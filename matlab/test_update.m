% test_update.m
% Test script for the update function

function test_update()
    disp('Running test_update...');

    % Test 1: Standard Velocity Verlet Update
    np = 2;
    nd = 2;
    pos = [0.0, 0.0; 0.0, 0.0];
    vel = [1.0, -1.0; 2.0, 0.0];
    f = [2.0, 0.0; 0.0, 4.0];
    acc = [0.0, 0.0; 0.0, 0.0];
    mass = 2.0;
    dt = 0.1;

    [new_pos, new_vel, new_acc] = update(np, nd, pos, vel, f, acc, mass, dt);

    expected_pos = [0.1, -0.1; 0.2, 0.0];
    expected_acc = [1.0, 0.0; 0.0, 2.0];
    expected_vel = [1.05, -1.0; 2.0, 0.1];

    assert(norm(new_pos - expected_pos) < 1e-6, 'Test 1 Failed: Position mismatch');
    assert(norm(new_vel - expected_vel) < 1e-6, 'Test 1 Failed: Velocity mismatch');
    assert(norm(new_acc - expected_acc) < 1e-6, 'Test 1 Failed: Acceleration mismatch');
    disp('Test 1 Passed.');

    % Test 2: Zero values
    np = 1;
    nd = 3;
    pos = [0.0; 0.0; 0.0];
    vel = [0.0; 0.0; 0.0];
    f = [0.0; 0.0; 0.0];
    acc = [0.0; 0.0; 0.0];
    mass = 1.0;
    dt = 0.5;

    [new_pos, new_vel, new_acc] = update(np, nd, pos, vel, f, acc, mass, dt);

    expected_pos = [0.0; 0.0; 0.0];
    expected_acc = [0.0; 0.0; 0.0];
    expected_vel = [0.0; 0.0; 0.0];

    assert(norm(new_pos - expected_pos) < 1e-6, 'Test 2 Failed: Position mismatch');
    assert(norm(new_vel - expected_vel) < 1e-6, 'Test 2 Failed: Velocity mismatch');
    assert(norm(new_acc - expected_acc) < 1e-6, 'Test 2 Failed: Acceleration mismatch');
    disp('Test 2 Passed.');

    % Test 3: Negative values and larger dt
    np = 1;
    nd = 2;
    pos = [1.0; -1.0];
    vel = [-0.5; 0.5];
    f = [-2.0; 2.0];
    acc = [0.5; -0.5];
    mass = 0.5;
    dt = 2.0;

    [new_pos, new_vel, new_acc] = update(np, nd, pos, vel, f, acc, mass, dt);

    expected_pos = [1.0; -1.0];
    expected_acc = [-4.0; 4.0];
    expected_vel = [-4.0; 4.0];

    assert(norm(new_pos - expected_pos) < 1e-6, 'Test 3 Failed: Position mismatch');
    assert(norm(new_vel - expected_vel) < 1e-6, 'Test 3 Failed: Velocity mismatch');
    assert(norm(new_acc - expected_acc) < 1e-6, 'Test 3 Failed: Acceleration mismatch');
    disp('Test 3 Passed.');

    disp('All tests passed successfully.');
end
