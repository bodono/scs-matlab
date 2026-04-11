classdef workspace < matlab.unittest.TestCase

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
            testCase.data.c = randn(n,1);
            testCase.cones.l = m;
            % Construct feasible b = A*x_feas + s_feas, s_feas >= 0
            x_feas = randn(n,1);
            s_feas = ones(m,1);
            testCase.data.b = testCase.data.A * x_feas + s_feas;
        end
    end

    methods (Test)
        function test_workspace_reuse(testCase, solver)
            pars = workspace.solver_pars(solver);
            pars.verbose = 0;

            % One-shot solve for reference
            [x_ref,~,~,info_ref] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info_ref.status, 'solved')

            % Workspace solve (same problem, should match)
            work = scs_init(testCase.data, testCase.cones, pars);
            [x_ws,~,~,info_ws] = scs_solve(work);
            testCase.verifyEqual(info_ws.status, 'solved')
            testCase.verifyEqual(x_ws, x_ref, 'RelTol', 1e-6)

            scs_finish(work);
        end

        function test_update_b(testCase, solver)
            pars = workspace.solver_pars(solver);
            pars.verbose = 0;

            work = scs_init(testCase.data, testCase.cones, pars);
            [~,~,~,info1] = scs_solve(work);
            testCase.verifyEqual(info1.status, 'solved')

            % Update b with a new feasible b and solve again
            rng(5678)
            x_feas = randn(size(testCase.data.c));
            s_feas = ones(size(testCase.data.b));
            b_new = testCase.data.A * x_feas + s_feas;
            scs_update(work, b_new, []);
            [~,~,~,info2] = scs_solve(work);
            testCase.verifyEqual(info2.status, 'solved')

            scs_finish(work);
        end

        function test_update_c(testCase, solver)
            pars = workspace.solver_pars(solver);
            pars.verbose = 0;

            work = scs_init(testCase.data, testCase.cones, pars);
            [~,~,~,info1] = scs_solve(work);
            testCase.verifyEqual(info1.status, 'solved')

            % Update c and solve again
            rng(5678)
            c_new = randn(size(testCase.data.c));
            scs_update(work, [], c_new);
            [~,~,~,info2] = scs_solve(work);
            testCase.verifyEqual(info2.status, 'solved')

            scs_finish(work);
        end

        function test_warm_start(testCase, solver)
            pars = workspace.solver_pars(solver);
            pars.verbose = 0;

            work = scs_init(testCase.data, testCase.cones, pars);
            [x,y,s,info1] = scs_solve(work);
            testCase.verifyEqual(info1.status, 'solved')

            % Re-solve with warm-start (same problem, should converge fast)
            warm.x = x;
            warm.y = y;
            warm.s = s;
            [~,~,~,info2] = scs_solve(work, warm);
            testCase.verifyEqual(info2.status, 'solved')
            testCase.verifyLessThanOrEqual(info2.iter, 25)

            scs_finish(work);
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
