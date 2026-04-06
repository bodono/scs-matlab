classdef mixed_cones < matlab.unittest.TestCase

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
            nz = 2;  % zero cone
            nl = 3;  % LP cone
            nq = 4;  % SOC of size 4
            ns_dim = 2;  % 2x2 SDP
            ns = ns_dim * (ns_dim + 1) / 2;
            m = nz + nl + nq + ns;
            testCase.data.A = sparse(randn(m, n));
            testCase.data.b = randn(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.z = nz;
            testCase.cones.l = nl;
            testCase.cones.q = nq;
            testCase.cones.s = ns_dim;
        end
    end

    methods (Test)
        function test_mixed(testCase, use_indirect)
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
