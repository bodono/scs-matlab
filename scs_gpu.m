function [x, y, s, info] = scs_gpu(data, cone, params)
% Operator-splitting method for solving cone problems (GPU indirect)
%
% This implements a cone solver. It solves:
%
% min. 0.5 * x'Px + c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% This uses the GPU indirect linear equation solver version of SCS.
%
% To use, set params.gpu = true when calling scs().
%
error ('scs_gpu mexFunction not found') ;
