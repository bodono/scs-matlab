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
            testCase.data.c = randn(n, 1);
            % Construct b so the problem is provably feasible:
            % pick x_feas, set s_feas in cone interior, b = A*x + s
            x_feas = randn(n, 1);
            s_feas = zeros(m, 1);
            % zero cone: s = 0
            % LP cone: s > 0
            s_feas(nz+1:nz+nl) = ones(nl, 1);
            % SOC: s_0 > ||s_rest||
            s_feas(nz+nl+1) = 2;
            s_feas(nz+nl+2:nz+nl+nq) = 0.5 * ones(nq-1, 1);
            % SDP (2x2): s = svec(I) = [1; 0; 1]
            s_feas(nz+nl+nq+1:m) = [1; 0; 1];
            testCase.data.b = testCase.data.A * x_feas + s_feas;
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
