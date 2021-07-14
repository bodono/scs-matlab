classdef basic < matlab.unittest.TestCase

    properties
        data
        cones
        %solver = {@scs_direct, @scs_indirect}
    end
    
    properties (TestParameter)
        solver = {@scs_direct, @scs_indirect}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Create Problem
            m = 9;
            n = 3;
            randn('seed',9)
            testCase.data.A = sparse(randn(m,n));
            testCase.data.b = randn(m,1);
            testCase.data.c = randn(n,1);
            testCase.cones.l = m;
        end
    end


    methods (Test)
        function test_random(testCase, solver)
            [x,y,s,info] = solver(testCase.data,testCase.cones,[]);
            testCase.verifyEqual(info.status, 'Solved')
            % warm-start test
            testCase.data.x = x;
            testCase.data.y = y;
            testCase.data.s = s;
            [~,~,~,info] = scs_direct(testCase.data,testCase.cones,[]);
            testCase.verifyEqual(info.status, 'Solved BAD')
            testCase.verifyEqual(info.iter, 0)
        end
   end
end

