classdef string_params < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        use_indirect = {false, true}
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
        function test_random(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.acceleration_lookback = 10;
            pars.write_data_filename = sprintf('data_dump_indirect_%d', use_indirect);
            %pars.log_csv_filename = sprintf('log_indirect_%d.csv', use_indirect);
            % test other params added
            pars.time_limit_secs = 11.0;
            pars.adaptive_scale = false;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')

            pars.time_limit_secs = 0.000001;
            pars.adaptive_scale = true;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved (inaccurate - reached time_limit_secs)')
        end
    end
end

