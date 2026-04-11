classdef info_fields < matlab.unittest.TestCase

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'matlab_ldl', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            rng(1234)
            m = 9;
            n = 3;
            testCase.data.A = sparse(randn(m,n));
            testCase.data.b = randn(m,1);
            testCase.data.c = randn(n,1);
            testCase.cones.l = m;
        end
    end

    methods (Test)
        function test_all_info_fields_present(testCase, solver)
            pars = info_fields.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);

            expected_fields = {'iter', 'status', 'pobj', 'dobj', ...
                'res_pri', 'res_dual', 'res_infeas', 'res_unbdd_a', ...
                'scale', 'status_val', 'res_unbdd_p', 'gap', ...
                'setup_time', 'solve_time', 'scale_updates', 'comp_slack'};
            for i = 1:length(expected_fields)
                testCase.verifyTrue(isfield(info, expected_fields{i}), ...
                    sprintf('Missing info field: %s', expected_fields{i}))
            end
        end

        function test_info_types(testCase, solver)
            pars = info_fields.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);

            testCase.verifyTrue(ischar(info.status))
            testCase.verifyTrue(isnumeric(info.iter))
            testCase.verifyTrue(isnumeric(info.pobj))
            testCase.verifyTrue(isnumeric(info.dobj))
            testCase.verifyTrue(isnumeric(info.setup_time))
            testCase.verifyTrue(isnumeric(info.solve_time))
            testCase.verifyGreaterThan(info.setup_time, 0)
            testCase.verifyGreaterThan(info.solve_time, 0)
        end

        function test_solved_status_val(testCase, solver)
            pars = info_fields.solver_pars(solver);
            pars.verbose = 0;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(info.status_val, 1)
        end

        function test_output_dimensions(testCase, solver)
            pars = info_fields.solver_pars(solver);
            pars.verbose = 0;
            [x,y,s,~] = scs(testCase.data,testCase.cones,pars);
            [m, n] = size(testCase.data.A);
            testCase.verifySize(x, [n, 1])
            testCase.verifySize(y, [m, 1])
            testCase.verifySize(s, [m, 1])
        end

        function test_duality_gap(testCase, solver)
            pars = info_fields.solver_pars(solver);
            pars.verbose = 0;
            pars.eps_abs = 1e-9;
            pars.eps_rel = 1e-9;
            [~,~,~,info] = scs(testCase.data,testCase.cones,pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyLessThan(abs(info.pobj - info.dobj), 1e-4)
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
