classdef infeasible < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        use_indirect = {false, true}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Infeasible LP: min 0 s.t. x = 1, x = -1
            testCase.data.A = sparse([1; -1]);
            testCase.data.b = [1; 1];
            testCase.data.c = [0];
            testCase.cones.z = 2;
        end
    end

    methods (Test)
        function test_infeasible(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'infeasible'))
        end
    end
end
