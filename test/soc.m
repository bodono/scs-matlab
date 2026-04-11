classdef soc < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'direct', 'indirect', 'matlab_ldl'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(1234)
            n = 5;
            m = n + 1;
            testCase.data.A = sparse([-ones(1,n); randn(n,n)]);
            testCase.data.b = zeros(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.q = m;
        end
    end

    methods (Test)
        function test_soc(testCase, solver)
            pars = soc.solver_pars(solver);
            pars.verbose = 0;
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
