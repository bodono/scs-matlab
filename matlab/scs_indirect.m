function [x, y, s, info] = scs_indirect(data, cone, params)
% Operator-splitting method for solving cone problems (indirect)
%
% This implements a cone solver. It solves:
%
% min. 0.5 * x'Px + c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% This uses the indirect (conjugate gradient) linear equation solver
% version of SCS.
%
% K is product of cones in this particular order:
% zero cone, lp cone, box cone, second order cone(s), semi-definite
% cone(s), complex semi-definite cone(s), primal exponential cones,
% dual exponential cones, power cones
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above.
% data.P is optional (set to [] or omit for LP/SOCP/SDP).
%
% cone struct must consist of:
% cone.z, length of zero cone (for equality constraints)
% cone.l, length of lp cone
% cone.bl, cone.bu, lower and upper bounds for box cone
% cone.q, array of SOC lengths
% cone.s, array of SD lengths
% cone.cs, array of complex SD lengths
% cone.ep, number of primal exp cones
% cone.ed, number of dual exp cones
% cone.p, array of power cone parameters
%
% Optional fields in the params struct are:
%   alpha                  : Douglas-Rachford relaxation parameter, between (0,2)
%   rho_x                  : primal constraint scaling factor
%   max_iters              : maximum number of iterations
%   eps_abs                : absolute convergence tolerance
%   eps_rel                : relative convergence tolerance
%   eps_infeas             : infeasibility tolerance
%   verbose                : verbosity level (0 or 1)
%   normalize              : heuristic data rescaling (0 or 1)
%   scale                  : initial dual scaling factor
%   adaptive_scale         : whether to adaptively update scale (0 or 1)
%   acceleration_lookback  : memory for Anderson acceleration (0 to disable)
%   acceleration_interval  : interval to apply acceleration
%   time_limit_secs        : time limit in seconds
%   write_data_filename    : if set, dump raw problem data to file
%   log_csv_filename       : if set, log progress to csv file
%
% to warm-start the solver add guesses for (x, y, s) to the data struct
%
error ('scs_indirect mexFunction not found') ;
