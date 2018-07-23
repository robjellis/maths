function vec_perm = permute_vec(vec)

% create a random permutation of vec

vec = vec(:);

nv = numel(vec);

% index
ind = randperm(nv);

% permute the vector
vec_perm = vec(ind);