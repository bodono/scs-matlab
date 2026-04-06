classdef qp < matlab.unittest.TestCase

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
            testCase.data.b = randn(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.l = m;
        end
    end

    methods (Test)
        function test_qp(testCase, use_indirect)
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

        function test_qp_lower_triangular_P(testCase, use_indirect)
            % scs.m should handle lower triangular P via symmetrization
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            % Pass only the lower triangle of the symmetric P
            P_lower = tril(testCase.data.P);
            testCase.data.P = sparse(P_lower);
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
        end
    end
end
