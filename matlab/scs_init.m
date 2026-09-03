function work = scs_init(data, K, pars)
% SCS_INIT  Initialize SCS workspace for repeated solves.
%
%   work = scs_init(data, K, pars)
%
%   Factorizes A (and P if present) and returns a workspace handle.
%   Use scs_solve(work) to solve, scs_update(work, b, c) to change
%   b/c without re-factorizing, and scs_finish(work) to free memory.
%
%   This is useful when solving a sequence of problems where only b
%   and/or c change, avoiding repeated factorization.
%
%   See also: scs_solve, scs_update, scs_finish, scs

if nargin < 3
    pars = [];
end

data = scs_prepare_data(data);

work.backend = scs_backend(pars);


work.id = feval(work.backend, 'init', data, K, pars);
