% Find the K channels with the smallest value Z
%
% [A, A_c] = find_good_channels(Z, K)
% 
% INPUT
%   Z       1 x BLOCKLENGTH vector
%   K       scalar
%
% OUTPUT
%   A       1 x BLOCKLENGTH logical vector
%   A_c     1 x BLOCKLENGTH logical vector
function [A, A_c] = find_good_channels(Z, K)
BLOCKLENGTH = length(Z);

[~, sorted_indices] = sort(Z);

A = false(1, BLOCKLENGTH);
A(sorted_indices(1:K)) = true;

A_c = false(1, BLOCKLENGTH);
A_c(sorted_indices(K+1:end)) = true;

end
