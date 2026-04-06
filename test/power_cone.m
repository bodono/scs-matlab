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
            % Simple power cone problem: min x3 s.t. cone constraints
            % With identity A, cone constraints directly bound x.
            n = 3;
            m = 6; % 2 power cone triples
            testCase.data.A = sparse([eye(3); eye(3)]);
            testCase.data.c = [0; 0; 1];
            % s = b - x, so cone constraints bound x directly
            % Primal power cone p=0.3: s1^0.3 * s2^0.7 >= |s3|
            % Dual power cone p=-0.7: dual cone constraint
            testCase.data.b = [1; 1; 0; 1; 1; 0];
            testCase.cones.p = [0.3, -0.7];
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
