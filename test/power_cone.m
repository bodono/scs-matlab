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
            % min x3 s.t. x >= 0 (LP), power cone constraint
            % LP rows: -x + s_lp = 0, s_lp >= 0 => x >= 0
            % Power rows: x + s_pow = [1;1;0]
            %   s_pow = [1-x1; 1-x2; -x3]
            %   (1-x1)^0.3 * (1-x2)^0.7 >= |x3|
            % With x >= 0: 0 <= x1 <= 1, 0 <= x2 <= 1, |x3| <= 1
            % Optimal: x3 = 0, obj = 0 (x3 >= 0 from LP)
            testCase.data.A = sparse([-eye(3); eye(3)]);
            testCase.data.b = [0; 0; 0; 1; 1; 0];
            testCase.data.c = [0; 0; 1];
            testCase.cones.l = 3;
            testCase.cones.p = 0.3;
        end
    end

    methods (Test)
        function test_power(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end
    end
end
