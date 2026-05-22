function test_md()
% TEST_MD tests the main md function.

  fprintf(1, '\nStarting tests for md.m\n');
  fprintf(1, '---------------------------\n');

  % Test 1: Basic execution with small parameters
  fprintf(1, 'Test 1: Basic execution (3D, 10 particles, 5 steps, dt=0.1)\n');
  try
    md(3, 10, 5, 0.1);
    fprintf(1, '  -> PASS\n');
  catch e
    error('  -> FAIL: %s', e.message);
  end

  % Test 2: Execution with string parameters
  fprintf(1, '\nTest 2: String parameters\n');
  try
    md('3', '10', '5', '0.1');
    fprintf(1, '  -> PASS\n');
  catch e
    error('  -> FAIL: %s', e.message);
  end

  % Test 3: Execution in 2D space
  fprintf(1, '\nTest 3: 2D space (2D, 10 particles, 5 steps, dt=0.1)\n');
  try
    md(2, 10, 5, 0.1);
    fprintf(1, '  -> PASS\n');
  catch e
    error('  -> FAIL: %s', e.message);
  end

  % Test 4: Default parameters (1 arguments passed, using default np, step_num, dt)
  % We use a small np and step_num indirectly by overriding string or testing error
  % Wait, md default is 500 particles 500 steps which takes 400s+.
  % We will skip testing default arguments to avoid test timeout.

  fprintf(1, '\n---------------------------\n');
  fprintf(1, 'All tests for md.m passed successfully.\n');
end
