classdef power_cone < matlab.unittest.TestCase

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
            n = 4;
            m = 6; % 2 power cone triples
            testCase.data.A = sparse(randn(m, n));
            testCase.data.c = randn(n, 1);
            testCase.cones.p = [0.3, -0.7]; % 1 primal, 1 dual
            % Construct feasible b = A*x_feas + s_feas
            x_feas = randn(n, 1);
            s_feas = zeros(m, 1);
            % Primal power cone p=0.3: x1^0.3 * x2^0.7 >= |x3|
            s_feas(1:3) = [1; 1; 0.5];
            % Dual power cone p=-0.7 (alpha=0.7):
            s_feas(4:6) = [0.7; 0.3; 0];
            testCase.data.b = testCase.data.A * x_feas + s_feas;
        end
    end

    methods (Test)
        function test_power(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'solved'))
        end
    end
end
