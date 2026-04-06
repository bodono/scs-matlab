function [x, y, s, info] = scs_solve(work, warm)
% SCS_SOLVE  Solve using an initialized SCS workspace.
%
%   [x, y, s, info] = scs_solve(work)
%   [x, y, s, info] = scs_solve(work, warm)
%
%   Solves the problem using the workspace from scs_init. Optionally
%   pass a warm-start struct with fields x, y, s from a previous solve.
%
%   See also: scs_init, scs_update, scs_finish

if nargin < 2
    [x, y, s, info] = feval(work.backend, 'solve');
else
    [x, y, s, info] = feval(work.backend, 'solve', warm);
end
