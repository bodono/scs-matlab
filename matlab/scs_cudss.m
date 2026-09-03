function varargout = scs_cudss(varargin)
% SCS_CUDSS  cuDSS GPU direct-solver backend (source build only).
%
%   Solves  min. 0.5 x'Px + c'x  s.t.  Ax + s = b,  s in K  on the GPU
%   through the NVIDIA cuDSS direct linear solver. It is not part of the
%   precompiled toolbox: build from source with SCS_BUILD_GPU=true (requires
%   the CUDA toolkit and cuDSS; see README) and the MEX takes precedence
%   over this stub. Select it with pars.gpu = true when calling scs.

error(['scs_cudss MEX not found: build from source with SCS_BUILD_GPU=true ' ...
       '(requires CUDA + cuDSS)']);
