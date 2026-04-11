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

if isfield(data, 'P')
    data.P = sparse(data.P);
    if (~istriu(data.P))
        data.P = triu(data.P + data.P') / 2;
    end
end

if isfield(data, 'A')
    data.A = sparse(data.A);
end

if size(data.b, 2) > 1
    data.b = data.b(:);
end

if size(data.c, 2) > 1
    data.c = data.c(:);
end

assert(size(data.A, 1) == size(data.b, 1), "A and b shape mismatch")
assert(size(data.A, 2) == size(data.c, 1), "A and c shape mismatch")

if isfield(data, 'P')
    assert(size(data.P, 1) == size(data.P, 2), "P is not square")
    assert(size(data.P, 1) == size(data.c, 1), "P and c shape mismatch")
end

if (isfield(pars,'use_indirect') && pars.use_indirect)
    work.backend = 'scs_indirect';
elseif (isfield(pars,'gpu') && pars.gpu)
    work.backend = 'scs_gpu';
elseif (isfield(pars,'dense') && pars.dense)
    work.backend = 'scs_dense';
elseif (isfield(pars,'use_qdldl') && pars.use_qdldl)
    work.backend = 'scs_direct';
else
    work.backend = 'scs_matlab_direct';
end

work.n = size(data.A, 2);
work.m = size(data.A, 1);

feval(work.backend, 'init', data, K, pars);
