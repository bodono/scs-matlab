classdef update_validation < matlab.unittest.TestCase

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
            testCase.data.A = sparse(randn(m, n));
            testCase.data.b = randn(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.l = m;
        end
    end

    methods (Test)
        function test_rejects_wrong_b_length(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            testCase.verifyError( ...
                @() scs_update(work, testCase.data.b(1:end-1), []), ...
                'scs:updateInvalidB')
        end

        function test_rejects_wrong_c_length(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            testCase.verifyError( ...
                @() scs_update(work, [], [testCase.data.c; 0]), ...
                'scs:updateInvalidC')
        end

        function test_rejects_sparse_update(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            testCase.verifyError( ...
                @() scs_update(work, sparse(testCase.data.b), []), ...
                'scs:updateInvalidB')
        end

        function test_rejects_non_double_update(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            testCase.verifyError( ...
                @() scs_update(work, [], single(testCase.data.c)), ...
                'scs:updateInvalidC')
        end

        function test_rejects_complex_update(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            testCase.verifyError( ...
                @() scs_update(work, testCase.data.b + 1i, []), ...
                'scs:updateInvalidB')
        end

        function test_accepts_row_vectors(testCase, solver)
            work = testCase.init_workspace(solver);
            cleanup = onCleanup(@() scs_finish(work)); %#ok<NASGU>

            scs_update(work, testCase.data.b', testCase.data.c')
            [~, ~, ~, info] = scs_solve(work);
            testCase.verifyEqual(info.status, 'solved')
        end
    end

    methods
        function work = init_workspace(testCase, solver)
            pars = struct('verbose', 0);
            if strcmp(solver, 'qdldl'), pars.use_qdldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
            work = scs_init(testCase.data, testCase.cones, pars);
        end
    end
end
