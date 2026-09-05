function varargout = scs(data, K, pars)
% scs 3.3.0
% for version call: scs_version()

nargoutchk(0, 4);

if nargin < 3
    pars = [];
end

data = scs_prepare_data(data);

[varargout{1:nargout}] = feval(scs_backend(pars), data, K, pars);
