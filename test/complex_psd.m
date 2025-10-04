classdef complex_psd < matlab.unittest.TestCase
    
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
            testCase.cones.cs = sqrt(m);
        end
    end

    methods (Test)
        function test_random(testCase, use_indirect)
            pars.use_indirect = use_indirect;
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

