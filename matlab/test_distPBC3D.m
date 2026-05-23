classdef test_distPBC3D < matlab.unittest.TestCase

    methods (Test)

        function testInsideBox(testCase)
            L = 10.0;
            vec = [1.0, 2.0, 3.0];
            expected = [1.0, 2.0, 3.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

        function testExceedsPositiveBoundary(testCase)
            L = 10.0;
            vec = [6.0, 2.0, 3.0];
            expected = [-4.0, 2.0, 3.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);

            vec = [1.0, 7.0, 3.0];
            expected = [1.0, -3.0, 3.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);

            vec = [1.0, 2.0, 8.0];
            expected = [1.0, 2.0, -2.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

        function testExceedsNegativeBoundary(testCase)
            L = 10.0;
            vec = [-6.0, -2.0, -3.0];
            expected = [4.0, -2.0, -3.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);

            vec = [-1.0, -7.0, -3.0];
            expected = [-1.0, 3.0, -3.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);

            vec = [-1.0, -2.0, -8.0];
            expected = [-1.0, -2.0, 2.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

        function testExactlyOnBoundary(testCase)
            L = 10.0;
            vec = [5.0, -5.0, 5.0];
            expected = [5.0, -5.0, 5.0];
            actual = distPBC3D(vec, L);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

    end
end