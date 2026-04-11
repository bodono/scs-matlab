classdef matlab_direct < matlab.unittest.TestCase

    properties
        data
        cones
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
        function test_basic_lp(testCase)
            pars.matlab_ldl = true;
            pars.verbose = 0;
            [x,y,s,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyTrue(contains(info.lin_sys_solver, 'matlab'))

            % warm-start test
            testCase.data.x = x;
            testCase.data.y = y;
            testCase.data.s = s;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyLessThanOrEqual(info.iter, 25)
        end

        function test_qp(testCase)
            rng(5678)
            n = size(testCase.data.c, 1);
            P = randn(n, n);
            P = P * P';
            testCase.data.P = sparse(P);
            pars.matlab_ldl = true;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_soc(testCase)
            rng(1234)
            n = 5;
            m = n + 1;
            testCase.data.A = sparse([-ones(1,n); randn(n,n)]);
            testCase.data.b = zeros(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones = struct('q', m);
            pars.matlab_ldl = true;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_matches_direct(testCase)
            % Verify MATLAB LDL produces the same solution as QDLDL
            pars_qdldl.verbose = 0;
            [x_ref,~,~,info_ref] = scs(testCase.data,testCase.cones,pars_qdldl);
            testCase.verifyEqual(info_ref.status, 'solved')

            pars_matlab.matlab_ldl = true;
            pars_matlab.verbose = 0;
            [x_ml,~,~,info_ml] = scs(testCase.data,testCase.cones,pars_matlab);
            testCase.verifyEqual(info_ml.status, 'solved')
            testCase.verifyEqual(x_ml, x_ref, 'RelTol', 1e-4)
        end

        function test_workspace(testCase)
            pars.matlab_ldl = true;
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

        function test_workspace_update(testCase)
            pars.matlab_ldl = true;
            pars.verbose = 0;

            % Construct a feasible problem
            rng(1234)
            m = 9; n = 3;
            testCase.data.A = sparse(randn(m,n));
            x_feas = randn(n,1);
            s_feas = ones(m,1);
            testCase.data.b = testCase.data.A * x_feas + s_feas;
            testCase.data.c = randn(n,1);
            testCase.cones = struct('l', m);

            work = scs_init(testCase.data, testCase.cones, pars);
            [~,~,~,info1] = scs_solve(work);
            testCase.verifyEqual(info1.status, 'solved')

            % Update b and re-solve
            rng(5678)
            x_feas2 = randn(n,1);
            b_new = testCase.data.A * x_feas2 + s_feas;
            scs_update(work, b_new, []);
            [~,~,~,info2] = scs_solve(work);
            testCase.verifyEqual(info2.status, 'solved')

            scs_finish(work);
        end

        function test_warm_start(testCase)
            pars.matlab_ldl = true;
            pars.verbose = 0;

            % Construct a feasible problem
            rng(1234)
            m = 9; n = 3;
            testCase.data.A = sparse(randn(m,n));
            x_feas = randn(n,1);
            s_feas = ones(m,1);
            testCase.data.b = testCase.data.A * x_feas + s_feas;
            testCase.data.c = randn(n,1);
            testCase.cones = struct('l', m);

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
end
