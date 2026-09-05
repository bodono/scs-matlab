function varargout = scs_solve(work, warm)
% SCS_SOLVE  Solve using an initialized SCS workspace.
%
%   [x, y, s, info] = scs_solve(work)
%   [x, y, s, info] = scs_solve(work, warm)
%
%   Solves the problem using the workspace from scs_init. Optionally
%   pass a warm-start struct with fields x, y, s from a previous solve.
%
%   See also: scs_init, scs_update, scs_finish

nargoutchk(0, 4);

if nargin < 2
    [varargout{1:nargout}] = feval(work.backend, 'solve', work.id);
else
    [varargout{1:nargout}] = feval(work.backend, 'solve', work.id, warm);
end
