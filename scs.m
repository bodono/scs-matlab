function [ x, y, s, info ] = scs( varargin )
% scs 3.0.0
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

if (isfield(pars,'use_indirect') && pars.use_indirect)
    [  x, y, s, info  ] = scs_indirect( data, K, pars );
elseif (isfield(pars,'gpu') && pars.gpu)
    [  x, y, s, info  ] = scs_gpu( data, K, pars );
else
    [  x, y, s, info  ] = scs_direct( data, K, pars );
end
