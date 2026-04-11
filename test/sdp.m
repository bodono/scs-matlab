classdef sdp < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'qdldl', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(1234)
            sd = 3;
            m = sd * (sd + 1) / 2; % = 6
            n = 4;
            testCase.data.A = sparse(randn(m, n));
            testCase.data.c = randn(n, 1);
            testCase.cones.s = sd;
            % Construct feasible b = A*x_feas + s_feas
            x_feas = randn(n, 1);
            % s_feas = svec(I_3): scaled lower triangle of identity
            % [1, 0, 0; 0, 1, 0; 0, 0, 1] -> [1; 0; 0; 1; 0; 1]
            s_feas = [1; 0; 0; 1; 0; 1];
            testCase.data.b = testCase.data.A * x_feas + s_feas;
        end
    end

    methods (Test)
        function test_sdp(testCase, solver)
            pars = sdp.solver_pars(solver);
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

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'qdldl'), pars.use_qdldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end
    end
end
