classdef power_cone < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        use_indirect = {false, true}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % min x3 s.t. (1-x1)^0.3 * (1-x2)^0.7 >= |x3|, x1<=1, x2<=1
            % Feasible: x=[0;0;0], optimal: x=[0;0;-1]
            testCase.data.A = speye(3);
            testCase.data.b = [1; 1; 0];
            testCase.data.c = [0; 0; 1];
            testCase.cones.p = 0.3;
        end
    end

    methods (Test)
        function test_power(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'solved'))
        end
    end
end
