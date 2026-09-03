classdef workspace_handles < matlab.unittest.TestCase

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
    end

    methods (Test)
        function test_same_backend_workspaces_are_independent(testCase, solver)
            data1.A = sparse(-1);
            data1.b = -1;
            data1.c = 1;
            cones1.l = 1;

            data2.A = -speye(2);
            data2.b = [-2; -3];
            data2.c = [1; 1];
            cones2.l = 2;

            pars = workspace_handles.solver_pars(solver);
            pars.verbose = 0;
            work1 = scs_init(data1, cones1, pars);
            cleanup1 = onCleanup(@() workspace_handles.finish_safely(work1)); %#ok<NASGU>
            work2 = scs_init(data2, cones2, pars);
            cleanup2 = onCleanup(@() workspace_handles.finish_safely(work2)); %#ok<NASGU>

            [x1, ~, ~, info1] = scs_solve(work1);
            [x2, ~, ~, info2] = scs_solve(work2);

            testCase.verifyEqual(info1.status, 'solved')
            testCase.verifyEqual(info2.status, 'solved')
            testCase.verifyEqual(x1, 1, 'AbsTol', 1e-3)
            testCase.verifyEqual(x2, [2; 3], 'AbsTol', 1e-3)
        end

        function test_finished_handle_cannot_target_new_workspace(testCase, solver)
            data.A = sparse(-1);
            data.b = -1;
            data.c = 1;
            cones.l = 1;
            pars = workspace_handles.solver_pars(solver);
            pars.verbose = 0;

            expired = scs_init(data, cones, pars);
            scs_finish(expired);
            current = scs_init(data, cones, pars);
            cleanup = onCleanup(@() workspace_handles.finish_safely(current)); %#ok<NASGU>

            testCase.verifyError(@() scs_solve(expired), 'scs:invalidWorkspace')
            [~, ~, ~, info] = scs_solve(current);
            testCase.verifyEqual(info.status, 'solved')
        end

        function test_live_handle_keeps_backend_loaded(testCase, solver)
            data.A = sparse(-1);
            data.b = -1;
            data.c = 1;
            cones.l = 1;
            pars = workspace_handles.solver_pars(solver);
            pars.verbose = 0;

            work = scs_init(data, cones, pars);
            cleanup = onCleanup(@() workspace_handles.finish_safely(work)); %#ok<NASGU>

            clear(work.backend)
            [x, ~, ~, info] = scs_solve(work);

            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(x, 1, 'AbsTol', 1e-3)
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'qdldl'), pars.use_qdldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end

        function finish_safely(work)
            try
                scs_finish(work);
            catch ME
                if ~strcmp(ME.identifier, 'scs:invalidWorkspace')
                    rethrow(ME)
                end
            end
        end
    end
end
