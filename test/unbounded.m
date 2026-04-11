classdef unbounded < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
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
        function test_unbounded(testCase, solver)
            pars = unbounded.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'unbounded'))
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'qdldl'), pars.use_qdldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end
    end
end
