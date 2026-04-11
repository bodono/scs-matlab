classdef matlab_ldl < matlab.unittest.TestCase
    % Cross-validation tests for the MATLAB LDL backend.
    % Compares solutions from default (matlab_ldl) against qdldl to ensure
    % the symmetrization, LDL extraction, and C-level solve are correct.

    methods (Test)
        function test_lp_cross_validate(testCase)
            % LP: no P matrix, simple linear cone
            rng(1234)
            m = 50; n = 20;
            data.A = sparse(randn(m, n));
            data.b = randn(m, 1);
            data.c = randn(n, 1);
            K.l = m;

            pars.verbose = 0;
            [x1, y1, ~, info1] = scs(data, K, pars);
            testCase.verifyEqual(info1.status, 'solved')

            pars.use_qdldl = true;
            [x2, y2, ~, info2] = scs(data, K, pars);
            testCase.verifyEqual(info2.status, 'solved')

            testCase.verifyEqual(x1, x2, 'RelTol', 1e-5, ...
                'matlab_ldl and qdldl should produce matching LP solutions')
            testCase.verifyEqual(y1, y2, 'RelTol', 1e-5)
        end

        function test_qp_cross_validate(testCase)
            % QP: with P matrix, exercises the P+R_x diagonal in KKT
            rng(2345)
            n = 15; m = 30;
            P = randn(n, n);
            P = sparse(P * P');
            data.A = sparse(randn(m, n));
            data.P = P;
            data.b = randn(m, 1);
            data.c = randn(n, 1);
            K.l = m;

            pars.verbose = 0;
            [x1, y1, ~, info1] = scs(data, K, pars);
            testCase.verifyEqual(info1.status, 'solved')

            pars.use_qdldl = true;
            [x2, y2, ~, info2] = scs(data, K, pars);
            testCase.verifyEqual(info2.status, 'solved')

            testCase.verifyEqual(x1, x2, 'RelTol', 1e-5, ...
                'matlab_ldl and qdldl should produce matching QP solutions')
            testCase.verifyEqual(y1, y2, 'RelTol', 1e-5)
        end

        function test_socp_cross_validate(testCase)
            % SOCP: second-order cone constraints
            rng(3456)
            n = 10;
            q_size = 15;
            data.A = sparse(randn(q_size, n));
            data.b = randn(q_size, 1);
            data.c = randn(n, 1);
            K.q = q_size;

            pars.verbose = 0;
            [x1, y1, ~, info1] = scs(data, K, pars);
            testCase.verifyEqual(info1.status, 'solved')

            pars.use_qdldl = true;
            [x2, y2, ~, info2] = scs(data, K, pars);
            testCase.verifyEqual(info2.status, 'solved')

            testCase.verifyEqual(x1, x2, 'RelTol', 1e-4, ...
                'matlab_ldl and qdldl should produce matching SOCP solutions')
            testCase.verifyEqual(y1, y2, 'RelTol', 1e-4)
        end

        function test_workspace_refactor_cross_validate(testCase)
            % Workspace with multiple updates triggers refactorization.
            % Verify both solvers produce the same results after refactor.
            rng(4567)
            n = 10; m = 20;
            data.A = sparse(randn(m, n));
            data.c = randn(n, 1);
            K.l = m;
            x_feas = randn(n, 1);
            s_feas = ones(m, 1);
            data.b = data.A * x_feas + s_feas;

            pars_default = struct('verbose', 0);
            pars_qdldl = struct('verbose', 0, 'use_qdldl', true);

            work_default = scs_init(data, K, pars_default);
            work_qdldl = scs_init(data, K, pars_qdldl);

            % Solve 3 times with different b vectors
            for i = 1:3
                rng(i * 100)
                x_feas = randn(n, 1);
                b_new = data.A * x_feas + s_feas;
                scs_update(work_default, b_new, []);
                scs_update(work_qdldl, b_new, []);

                [x1, ~, ~, info1] = scs_solve(work_default);
                [x2, ~, ~, info2] = scs_solve(work_qdldl);

                testCase.verifyEqual(info1.status, 'solved')
                testCase.verifyEqual(info2.status, 'solved')
                testCase.verifyEqual(x1, x2, 'RelTol', 1e-5, ...
                    sprintf('Iteration %d: workspace solutions should match', i))
            end

            scs_finish(work_default);
            scs_finish(work_qdldl);
        end

        function test_mixed_cones_cross_validate(testCase)
            % Mixed LP + SOC cones for a more complex KKT structure
            rng(5678)
            n = 10;
            m_lp = 5;
            m_soc = 8;
            m = m_lp + m_soc;
            data.A = sparse(randn(m, n));
            data.b = randn(m, 1);
            data.c = randn(n, 1);
            K.l = m_lp;
            K.q = m_soc;

            pars.verbose = 0;
            [x1, y1, ~, info1] = scs(data, K, pars);
            testCase.verifyEqual(info1.status, 'solved')

            pars.use_qdldl = true;
            [x2, y2, ~, info2] = scs(data, K, pars);
            testCase.verifyEqual(info2.status, 'solved')

            testCase.verifyEqual(x1, x2, 'RelTol', 1e-4, ...
                'matlab_ldl and qdldl should match on mixed cone problems')
            testCase.verifyEqual(y1, y2, 'RelTol', 1e-4)
        end
    end
end
