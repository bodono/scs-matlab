addpath ~/cvx
addpath ~/scs/matlab
addpath ~/ecos-matlab/bin

cvx_setup ~/cvx_license.dat

parameters.run_cvx = true;
parameters.save_results = true;
parameters.run_scs_direct = false;
parameters.run_scs_indirect = false;

solvers = {'mosek', 'gurobi', 'ecos', 'sdpt3', 'sedumi'}

for i = 1:length(solvers)
    parameters.cvx_use_solver = solvers{i}
    run_lasso_ex(parameters)
    run_portfolio_ex(parameters)
    run_rpca_ex(parameters)
    run_pnorm_ex(parameters)
    %run_l1logreg_ex(parameters)
end
