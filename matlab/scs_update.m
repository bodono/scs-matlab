function scs_update(work, b, c)
% SCS_UPDATE  Update b and/or c without re-factorizing.
%
%   scs_update(work, b_new, c_new)
%   scs_update(work, b_new, [])    % update only b
%   scs_update(work, [], c_new)    % update only c
%
%   Updates the b and/or c vectors in the workspace. Pass [] to leave
%   a vector unchanged. After updating, call scs_solve to re-solve.
%
%   See also: scs_init, scs_solve, scs_finish

if nargin < 3
    c = [];
end

feval(work.backend, 'update', work.id, b, c);
