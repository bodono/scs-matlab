classdef quad_box < matlab.unittest.TestCase
    
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
            P = randn(3, 3);
            P = P * P';
            testCase.data.A = sparse([zeros(1,3);randn(3,3)]);
            testCase.data.P = sparse(P);
            testCase.data.b = [1., 0., 0. 0.]';
            testCase.data.c = -ones(3,1);
            testCase.cones.bl = [0., 1., -2.];
            testCase.cones.bu = [1., 2., -1.];
        end
    end

    methods (Test)
        function test_random(testCase, use_indirect)
            pars.use_indirect = use_indirect;
            pars.acceleration_lookback = 10;
            [x,y,s,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            
            % warm-start test
            testCase.data.x = x;
            testCase.data.y = y;
            testCase.data.s = s;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(info.iter, 0)
        end
    end
end

