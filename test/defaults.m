classdef defaults < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
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
        function test_no_params(testCase)
            % Call scs with only 2 arguments (no params struct)
            [~,~,~,info] = scs(testCase.data,testCase.cones);
            testCase.verifyEqual(info.status, 'solved')
            % Verify that the default solver is MATLAB's LDL
            testCase.verifyEqual(info.lin_sys_solver, 'sparse-direct-matlab-ldl')
        end

        function test_empty_params(testCase)
            [~,~,~,info] = scs(testCase.data,testCase.cones,[]);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_workspace_no_params(testCase)
            % Call scs_init with only 2 arguments
            work = scs_init(testCase.data, testCase.cones);
            [~,~,~,info] = scs_solve(work);
            scs_finish(work);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_workspace_empty_params(testCase)
            % Call scs_init with empty params
            work = scs_init(testCase.data, testCase.cones, []);
            [~,~,~,info] = scs_solve(work);
            scs_finish(work);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_legacy_f_cone(testCase)
            % Verify that f and z cones are summed correctly
            data.A = sparse([1 0; 0 1; -1 0; 0 -1; 1 1]);
            data.b = [1; 1; 0; 0; 1];
            data.c = [-1; -1];
            % Define 2 zero cones via 'f' and 2 via 'z' (total 4)
            % plus one linear cone 'l' (total 5 rows)
            cones.f = 2;
            cones.z = 2;
            cones.l = 1;
            [~,~,~,info] = scs(data, cones, struct('verbose', 0));
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_dense_A_auto_sparsified(testCase, solver)
            % Pass dense A — scs.m should auto-convert to sparse
            testCase.data.A = full(testCase.data.A);
            pars = defaults.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_row_vector_b_c(testCase, solver)
            % Pass b and c as row vectors — scs.m should reshape
            testCase.data.b = testCase.data.b';
            testCase.data.c = testCase.data.c';
            pars = defaults.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
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
