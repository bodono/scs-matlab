function scs_finish(work)
% SCS_FINISH  Free SCS workspace memory.
%
%   scs_finish(work)
%
%   Frees the workspace allocated by scs_init. Must be called when
%   done to avoid memory leaks.
%
%   See also: scs_init, scs_solve, scs_update

feval(work.backend, 'finish');
