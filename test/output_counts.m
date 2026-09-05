classdef output_counts < matlab.unittest.TestCase
    methods (Test)
        function test_one_shot_short_outputs(testCase)
            data = struct('A', sparse(-1), 'b', -1, 'c', 1);
            cones = struct('l', 1);
            pars = struct('verbose', 0);

            % Exercise the MEX's zero-output frame and one-output copy
            % directly, as well as nargout forwarding by the public API.
            scs_direct(data, cones, pars);
            x = scs_direct(data, cones, pars);
            testCase.verifyEqual(x, 1, 'AbsTol', 1e-3);
            scs(data, cones, pars);
            x = scs(data, cones, pars);
            testCase.verifyEqual(x, 1, 'AbsTol', 1e-3);
        end
    end
end
