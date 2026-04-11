classdef string_params < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Create Problem
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
        function test_random(testCase, solver)
            pars = string_params.solver_pars(solver);
            pars.acceleration_lookback = 10;
            % test string parameters
            pars.write_data_filename = sprintf('data_dump_%s', solver);
            pars.log_csv_filename = sprintf('log_%s.csv', solver);
            % test other parameters added at the same time
            pars.time_limit_secs = 11.0;
            pars.adaptive_scale = true;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')

            pars.time_limit_secs = 0.000001;
            pars.adaptive_scale = false;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved (inaccurate - reached time_limit_secs)')
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
