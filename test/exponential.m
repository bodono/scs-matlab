classdef exponential < matlab.unittest.TestCase

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
            n = 4;
            m = 6; % 2 exp cone triples (1 primal + 1 dual)
            testCase.data.A = sparse(randn(m, n));
            testCase.data.b = randn(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.ep = 1;
            testCase.cones.ed = 1;
        end
    end

    methods (Test)
        function test_exponential(testCase, solver)
            pars = exponential.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'solved') || ...
                contains(info.status, 'infeasible'))
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
