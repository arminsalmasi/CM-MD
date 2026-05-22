function tests = test_initialize
  tests = functiontests(localfunctions);
end

function testCorrectDimensions(testCase)
  % Test to ensure initialize returns matrices with the correct dimensions
  np = 100;
  nd = 3;
  [pos, vel, acc] = initialize(np, nd);

  testCase.verifySize(pos, [nd, np]);
  testCase.verifySize(vel, [nd, np]);
  testCase.verifySize(acc, [nd, np]);
end

function testZeroInitialization(testCase)
  % Test that velocity and acceleration are initialized to zeros
  np = 50;
  nd = 2;
  [~, vel, acc] = initialize(np, nd);

  testCase.verifyEqual(vel, zeros(nd, np));
  testCase.verifyEqual(acc, zeros(nd, np));
end

function testPositionBounds(testCase)
  % Test that position is initialized between 0 and 10 (as per the code: 0.0, 10.0)
  np = 1000;
  nd = 3;
  [pos, ~, ~] = initialize(np, nd);

  testCase.verifyGreaterThanOrEqual(pos, 0.0);
  testCase.verifyLessThanOrEqual(pos, 10.0);
end

function testDeterministicOutput(testCase)
  % Test that calling initialize twice with same params yields same pos
  % due to hardcoded seed
  np = 10;
  nd = 3;
  [pos1, ~, ~] = initialize(np, nd);
  [pos2, ~, ~] = initialize(np, nd);

  testCase.verifyEqual(pos1, pos2);
end
