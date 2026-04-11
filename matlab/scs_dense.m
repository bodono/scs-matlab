function [x, y, s, info] = scs_dense(data, cone, params)
% Operator-splitting method for solving cone problems (dense direct)
%
% This implements a cone solver. It solves:
%
% min. 0.5 * x'Px + c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% This uses the dense direct (Gram matrix + Cholesky) linear solver.
% Best suited for problems where A is dense.
%
% To use, set params.dense = true when calling scs().
%
error ('scs_dense mexFunction not found') ;
