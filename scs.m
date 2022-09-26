function [ x, y, s, info ] = scs( varargin )
% scs 3.2.0
% for version call: scs_version()
data = varargin{1};
K = varargin{2};
if nargin >= 3
    pars = varargin{3};
else
    pars = [];
end

if isfield(data, 'P')
    data.P = sparse(data.P);
    if (~istriu(data.P))
        data.P = triu(data.P);
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
    [  x, y, s, info  ] = scs_indirect( data, K, pars );
elseif (isfield(pars,'gpu') && pars.gpu)
    [  x, y, s, info  ] = scs_gpu( data, K, pars );
else
    [  x, y, s, info  ] = scs_direct( data, K, pars );
end
