classdef sdp < matlab.unittest.TestCase

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
            sd = 3;
            m = sd * (sd + 1) / 2;
            n = 4;
            testCase.data.A = sparse(randn(m, n));
            testCase.data.c = randn(n, 1);
            % Construct feasible b = A*x_feas + s_feas
            x_feas = randn(n, 1);
            % s_feas = svec(I_3): scaled lower triangle of identity
            s_feas = [1; 0; 0; 1; 0; 1];
            testCase.data.b = testCase.data.A * x_feas + s_feas;
            testCase.cones.s = sd;
        end
    end

    methods (Test)
        function test_sdp(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [x,y,s,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')

            % warm-start test
            testCase.data.x = x;
            testCase.data.y = y;
            testCase.data.s = s;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyLessThanOrEqual(info.iter, 25)
        end
    end
end
