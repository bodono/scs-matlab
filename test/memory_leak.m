classdef memory_leak < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'matlab_ldl', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
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
        function test_repeated_scs(testCase, solver)
            % Solve many times in a loop. If persistent mxArrays leak,
            % MATLAB's memory usage will grow unboundedly.
            pars = memory_leak.solver_pars(solver);
            pars.verbose = 0;
            n_runs = 100;
            for i = 1:n_runs
                [~,~,~,info] = scs(testCase.data, testCase.cones, pars);
                testCase.verifyEqual(info.status, 'solved')
            end
        end

        function test_repeated_workspace(testCase, solver)
            % Workspace init/solve/finish cycle. Tests that scs_finish
            % properly frees all allocated memory.
            pars = memory_leak.solver_pars(solver);
            pars.verbose = 0;
            n_runs = 100;
            for i = 1:n_runs
                work = scs_init(testCase.data, testCase.cones, pars);
                [~,~,~,info] = scs_solve(work);
                testCase.verifyEqual(info.status, 'solved')
                scs_finish(work);
            end
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
