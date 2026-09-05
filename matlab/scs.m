function varargout = scs(data, K, pars)
% scs 3.3.0
% for version call: scs_version()

nargoutchk(0, 4);

if nargin < 3
    pars = [];
end

data = scs_prepare_data(data);

if isfield(pars, 'use_indirect') && pars.use_indirect
    [varargout{1:nargout}] = scs_indirect(data, K, pars);
elseif isfield(pars, 'gpu') && pars.gpu
    [varargout{1:nargout}] = scs_gpu(data, K, pars);
elseif isfield(pars, 'dense') && pars.dense
    [varargout{1:nargout}] = scs_dense(data, K, pars);
elseif isfield(pars, 'use_qdldl') && pars.use_qdldl
    [varargout{1:nargout}] = scs_direct(data, K, pars);
else
    [varargout{1:nargout}] = scs_matlab_direct(data, K, pars);
end
