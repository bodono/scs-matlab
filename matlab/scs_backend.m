function name = scs_backend(pars)
% SCS_BACKEND  Name of the MEX backend that a settings struct selects.
%
%   use_indirect  -> scs_indirect
%   gpu           -> scs_cudss (source build only; requires CUDA + cuDSS)
%   dense         -> scs_dense
%   use_qdldl     -> scs_direct
%   otherwise     -> scs_matlab_direct
%
%   Used by both scs and scs_init so the two can never disagree.

if isfield(pars, 'use_indirect') && pars.use_indirect
    name = 'scs_indirect';
elseif isfield(pars, 'gpu') && pars.gpu
    name = 'scs_cudss';
elseif isfield(pars, 'dense') && pars.dense
    name = 'scs_dense';
elseif isfield(pars, 'use_qdldl') && pars.use_qdldl
    name = 'scs_direct';
else
    name = 'scs_matlab_direct';
end
