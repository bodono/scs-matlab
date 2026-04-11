classdef infeasible < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'matlab_ldl', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Infeasible LP: min 0 s.t. x = 1, x = -1
            testCase.data.A = sparse([1; -1]);
            testCase.data.b = [1; 1];
            testCase.data.c = [0];
            testCase.cones.z = 2;
        end
    end

    methods (Test)
        function test_infeasible(testCase, solver)
            pars = infeasible.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyTrue(contains(info.status, 'infeasible'))
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'matlab_ldl'), pars.use_matlab_ldl = true; end
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end
    end
end
