function x = scs_matlab_ldl_solve(L, D, p, b)
%SCS_MATLAB_LDL_SOLVE Solve LDL-factored linear system.
%   x = scs_matlab_ldl_solve(L, D, p, b) solves A*x = b given the LDL
%   factorization [L, D, p] = ldl(A, 'vector'), i.e., A(p,p) = L*D*L'.
x = zeros(size(b));
x(p) = L' \ (D \ (L \ b(p)));
end
