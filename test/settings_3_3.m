classdef settings_3_3 < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
        regularization = {1e-10, -1e-10}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.data.A = sparse(-1);
            testCase.data.b = -1;
            testCase.data.c = 1;
            testCase.cones.l = 1;
        end
    end

    methods (Test)
        function test_new_settings_are_accepted(testCase, solver, regularization)
            pars = settings_3_3.solver_pars(solver);
            pars.verbose = 0;
            pars.adaptive_diag_scale = 0;
            pars.acceleration_type_1 = 0;
            pars.acceleration_regularization = regularization;
            pars.acceleration_relaxation = 0.9;

            [~, ~, ~, info] = scs(testCase.data, testCase.cones, pars);

            testCase.verifyEqual(info.status, 'solved')
        end

        function test_invalid_adaptive_diag_scale_reaches_core(testCase, solver)
            pars = settings_3_3.solver_pars(solver);
            pars.verbose = 0;
            pars.adaptive_diag_scale = 2;
            [~, ~, ~, info] = scs(testCase.data, testCase.cones, pars);
            testCase.verifyEqual(info.status, 'failure')
        end

        function test_invalid_acceleration_relaxation_reaches_core(testCase, solver)
            pars = settings_3_3.solver_pars(solver);
            pars.verbose = 0;
            pars.acceleration_relaxation = 3;
            [~, ~, ~, info] = scs(testCase.data, testCase.cones, pars);
            testCase.verifyEqual(info.status, 'failure')
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
