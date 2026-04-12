function [x, y, s, info] = scs(data, K, pars)
% scs 3.2.5
% for version call: scs_version()

if nargin < 3
    pars = [];
end

data = scs_prepare_data(data);

if isfield(pars, 'use_indirect') && pars.use_indirect
    [x, y, s, info] = scs_indirect(data, K, pars);
elseif isfield(pars, 'gpu') && pars.gpu
    [x, y, s, info] = scs_gpu(data, K, pars);
elseif isfield(pars, 'dense') && pars.dense
    [x, y, s, info] = scs_dense(data, K, pars);
elseif isfield(pars, 'use_qdldl') && pars.use_qdldl
    [x, y, s, info] = scs_direct(data, K, pars);
else
    [x, y, s, info] = scs_matlab_direct(data, K, pars);
end
