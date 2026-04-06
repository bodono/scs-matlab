classdef unbounded < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        use_indirect = {false, true}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Unbounded LP: min -x s.t. x >= 0
            % Ax + s = b, s >= 0: -x + s = 0 => s = x >= 0
            testCase.data.A = sparse([-1]);
            testCase.data.b = [0];
            testCase.data.c = [-1];
            testCase.cones.l = 1;
        end
    end

    methods (Test)
        function test_unbounded(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'unbounded'))
        end
    end
end
