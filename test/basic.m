classdef basic < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'direct', 'indirect', 'matlab_ldl'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Create Problem
            rng(1234)
            m = 9;
            n = 3;
            testCase.data.A = sparse(randn(m,n));
            testCase.data.b = randn(m,1);
            testCase.data.c = randn(n,1);
            testCase.cones.l = m;
        end
    end

    methods (Test)
        function test_random(testCase, solver)
            pars = basic.solver_pars(solver);
            pars.acceleration_lookback = 10;
            [x,y,s,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')

            % warm-start test
            testCase.data.x = x;
            testCase.data.y = y;
            testCase.data.s = s;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyLessThanOrEqual(info.iter, 25)
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
            if strcmp(solver, 'matlab_ldl'), pars.matlab_ldl = true; end
        end
    end
end
