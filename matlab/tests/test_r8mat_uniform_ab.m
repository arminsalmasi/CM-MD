function test_r8mat_uniform_ab()
  disp('Running tests for r8mat_uniform_ab...');
  addpath('..'); % Go up to the matlab folder to find the function

  % Test 1: Dimensions
  disp('Test 1: Check dimensions...');
  [r, seed] = r8mat_uniform_ab(2, 3, 0.0, 10.0, 12345);
  assert(all(size(r) == [2, 3]), 'Matrix dimensions are incorrect');

  % Test 2: Boundaries
  disp('Test 2: Check boundaries...');
  % Create a larger matrix to test boundaries better
  [r, seed] = r8mat_uniform_ab(100, 100, -5.0, 5.0, 12345);
  assert(all(r(:) >= -5.0), 'Values are below lower bound');
  assert(all(r(:) <= 5.0), 'Values are above upper bound');

  % Test 3: Seed update
  disp('Test 3: Check seed updates...');
  [r1, seed1] = r8mat_uniform_ab(2, 2, 0.0, 1.0, 12345);
  [r2, seed2] = r8mat_uniform_ab(2, 2, 0.0, 1.0, seed1);
  assert(seed1 ~= 12345, 'Seed was not updated');
  assert(seed2 ~= seed1, 'Seed was not updated on second call');
  assert(any(r1(:) ~= r2(:)), 'Successive matrices are identical');

  % Test 4: Reproducibility
  disp('Test 4: Check reproducibility...');
  [r3, seed3] = r8mat_uniform_ab(2, 2, 0.0, 1.0, 12345);
  assert(all(r1(:) == r3(:)), 'Matrices generated with same seed are different');
  assert(seed1 == seed3, 'Updated seeds are different for same initial seed');

  % Test 5: Error on seed 0
  disp('Test 5: Check error on seed=0...');
  try
    [~, ~] = r8mat_uniform_ab(2, 2, 0.0, 1.0, 0);
    error('Did not throw error on seed=0');
  catch e
    assert(strcmp(e.message, 'R8MAT_UNIFORM_AB - Fatal error!'), 'Unexpected error message: %s', e.message);
  end

  disp('All tests passed successfully!');
end

test_r8mat_uniform_ab();
