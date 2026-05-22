classdef test_timestamp < matlab.unittest.TestCase
    % Tests for timestamp function

    methods (Test)
        function testFormat(testCase)
            % Capture the output of timestamp function
            outputStr = evalc('timestamp()');

            % Remove any trailing whitespace/newline
            outputStr = strtrim(outputStr);

            % The output should match the datestr(..., 0) format:
            % dd-mmm-yyyy HH:MM:SS
            % E.g., '27-Dec-2014 12:34:56'

            % Check that output is not empty
            testCase.verifyNotEmpty(outputStr, 'Output should not be empty');

            % Check format using regex
            % Format: 2 digits, '-', 3 letters, '-', 4 digits, ' ', 2 digits, ':', 2 digits, ':', 2 digits
            expression = '^\d{2}-[A-Za-z]{3}-\d{4} \d{2}:\d{2}:\d{2}$';
            match = regexp(outputStr, expression, 'once');

            testCase.verifyNotEmpty(match, sprintf('Output "%s" does not match expected format', outputStr));
        end

        function testCurrentTime(testCase)
            % Check if the printed time is close to current time

            % Get current time
            t1 = now;

            % Capture output
            outputStr = evalc('timestamp()');
            outputStr = strtrim(outputStr);

            % Get current time right after
            t2 = now;

            % Try to parse the output time
            try
                parsedTime = datenum(outputStr);

                % Check if parsed time is between t1 and t2 (with some margin for formatting truncation)
                % datenum with 0 format gives second precision, so add/subtract a small tolerance
                tolerance = 2 / (24 * 60 * 60); % 2 seconds

                testCase.verifyTrue(parsedTime >= t1 - tolerance, 'Parsed time is earlier than expected');
                testCase.verifyTrue(parsedTime <= t2 + tolerance, 'Parsed time is later than expected');
            catch ME
                testCase.verifyFail(sprintf('Failed to parse output string as datenum: %s', ME.message));
            end
        end
    end
end
