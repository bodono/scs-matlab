classdef workspace_qp < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        use_indirect = {false, true}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(1234)
            n = 5;
            m = 10;
            P = randn(n, n);
            P = P * P';
            testCase.data.A = sparse(randn(m, n));
            testCase.data.P = sparse(P);
            testCase.data.c = randn(n, 1);
            testCase.cones.l = m;
            % Construct feasible b = A*x_feas + s_feas, s_feas >= 0
            x_feas = randn(n, 1);
            s_feas = ones(m, 1);
            testCase.data.b = testCase.data.A * x_feas + s_feas;
        end
    end

    methods (Test)
        function test_qp_workspace(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;

            % One-shot reference
            [x_ref,~,~,info_ref] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info_ref.status, 'solved')

            % Workspace solve
            work = scs_init(testCase.data, testCase.cones, pars);
            [x_ws,~,~,info_ws] = scs_solve(work);
            testCase.verifyEqual(info_ws.status, 'solved')
            testCase.verifyEqual(x_ws, x_ref, 'RelTol', 1e-6)

            scs_finish(work);
        end

        function test_update_both_b_and_c(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;

            work = scs_init(testCase.data, testCase.cones, pars);
            [~,~,~,info1] = scs_solve(work);
            testCase.verifyEqual(info1.status, 'solved')

            % Update both b and c with feasible data
            rng(5678)
            x_feas = randn(size(testCase.data.c));
            s_feas = ones(size(testCase.data.b));
            b_new = testCase.data.A * x_feas + s_feas;
            c_new = randn(size(testCase.data.c));
            scs_update(work, b_new, c_new);
            [~,~,~,info2] = scs_solve(work);
            testCase.verifyEqual(info2.status, 'solved')

            scs_finish(work);
        end

        function test_sequential_updates(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;

            work = scs_init(testCase.data, testCase.cones, pars);

            % Solve 5 times with different feasible b vectors
            for i = 1:5
                rng(i * 1000)
                x_feas = randn(size(testCase.data.c));
                s_feas = ones(size(testCase.data.b));
                b_new = testCase.data.A * x_feas + s_feas;
                scs_update(work, b_new, []);
                [~,~,~,info] = scs_solve(work);
                testCase.verifyEqual(info.status, 'solved')
            end

            scs_finish(work);
        end
    end
end
