function data = scs_prepare_data(data)
% SCS_PREPARE_DATA Normalize input matrices/vectors and validate dimensions.

if isfield(data, 'P')
    data.P = sparse(data.P);
    if ~istriu(data.P)
        data.P = triu(data.P + data.P') / 2;
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
