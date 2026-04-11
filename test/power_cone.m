classdef power_cone < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'matlab_ldl', 'indirect'}
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
        function test_power(testCase, solver)
            pars = power_cone.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'matlab_ldl'), pars.use_matlab_ldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end
    end
end
