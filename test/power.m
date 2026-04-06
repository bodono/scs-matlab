classdef power < matlab.unittest.TestCase

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
            testCase.data.b = randn(m, 1);
            testCase.data.c = randn(n, 1);
            testCase.cones.p = [0.3, -0.7]; % 1 primal, 1 dual
        end
    end

    methods (Test)
        function test_power(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'solved') || ...
                contains(info.status, 'infeasible'))
        end
    end
end
