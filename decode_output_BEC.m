% Decode the received vector in an efficient way (with dynamic programming).
%
% decoded_output = decode_output_BEC_naive(received_output, frozen_bits, A, A_c)
%
% INPUT
%   received_output     1 x BLOCKLENGTH vector
%   frozen_bits         1 x (BLOCKLENGTH - K) vector
%   A                   1 x BLOCKLENGTH logical vector
%   A_c                 1 x BLOCKLENGTH logical vector
%
% OUTPUT
%   decoded_output      1 x K vector
function decoded_output = decode_output_BEC(received_output, frozen_bits, A, A_c)

BLOCKLENGTH = length(received_output);
INITIAL_LEVEL = 1;
INITIAL_SHIFT_J = 0;

% Put frozen bits in a 1 x BLOCKLENGTH vector, at positions A_c, for easier
% access
frozen_bits_expanded = nan(1,BLOCKLENGTH);
frozen_bits_expanded(A_c) = frozen_bits;

% Decode bit by bit in an optimized way (reusing previous results)
LRs = nan(BLOCKLENGTH, log2(BLOCKLENGTH)+1);

decoded_output = nan(1, BLOCKLENGTH);
for j = 1:BLOCKLENGTH
    if(A_c(j))
        % If the bit is frozen, we don't need to compute anything
        decoded_output(j) = frozen_bits_expanded(j);
    else
        % To decode, first compute the likelihood ratio using the previously
        % decoded bits
        [LRs, L] = compute_lr(received_output, decoded_output(1:j-1), ...
            BLOCKLENGTH, j, LRs, INITIAL_SHIFT_J, INITIAL_LEVEL);
        
        % Then decide according to the lr
        decoded_output(j) = decide(L);
        
        % If we cannot recover (erasure), we stop decoding
        if(isnan(decoded_output(j)))
            break;
        end
    end
end

% Only return the information bits
decoded_output = decoded_output(A);

end


% Compute the likelihood ration using channels outputs and decoded bits.
% See Arikan formula 74 and 75
%
% We complete the table LRs to avoid recomputing the same likelihood ratio.
% 
% level indicates the recursion level (1 for decision level, log2(N)+1 for
% channel level).
% 
% shift_j allows to compute quickly where the local j is in the LRs table
% (abs_j = j + shift_j).
% 
% For convenience, return L := LRs(j)
%
% INPUT
%   y           1 x N vector (reduced by 2 at each recursive call)
%   u           1 x (j-1) vector (reduced by 2 at each recursive call)
%   N           scalar (reduced by 2 at each recursive call)
%   j           scalar (in a range {1, ..., N})
%   LRs         N x N_LEVELS array, with N_LEVELS := log2(N)+1
%   shift_j     scalar
%   level       scalar
%
% OUTPUT
%   LRs         N x N_LEVELS array, with N_LEVELS := log2(N)+1
%   L           scalar
function [LRs, L] = compute_lr(y, u, N, j, LRs, shift_j, level)

j_abs = shift_j + j;

% If the desired likelihood ratio already exists in the table, return
if(~isnan(LRs(j_abs, level)))
    L = LRs(j_abs, level);
    return;
end

% Recusion termination condition L(y) = W(y|0) / W(y|1)
% Note that only this part is specific to BEC
if(N == 1)
    if(y == 0)
        L = Inf;
    elseif(y == 1)
        L = 0;
    elseif(isnan(y)) % Erasure
        L = 1;
    else
        error('Unexpected value for y')
    end
    
    LRs(j_abs, level) = L;
    return;
end

% % Use formula 74 or 75 according to the parity of j
% First compute the recursive likelihood ratio (same are used)
if(mod(j,2)==1)
    j_even = j+1;
else
    j_even = j;
end

u_odd = u(1:2:j_even-2);
u_even = u(2:2:j_even-1);

[LRs, L1] = compute_lr(y(1:N/2), mod(u_odd + u_even, 2), N/2, j_even/2, ...
    LRs, shift_j, level+1);

[LRs, L2] = compute_lr(y(N/2+1:N), u_even, N/2, j_even/2, ...
    LRs, shift_j+N/2, level+1);

% Then combine the likelihood ratio
if(mod(j,2) == 1)
    % Use table decision for border cases 0/0, Inf/Inf (make diagram to understand)
    if((L1 == 0 && L2 == 0) || (isinf(L1) && isinf(L2)))
        LRs(j_abs, level) = Inf;
    elseif((L1 == 0 && isinf(L2)) || (isinf(L1) && L2 == 0))
        LRs(j_abs, level) = 0;
    elseif((L1 == 1 && isinf(L2)) || (isinf(L1) && L2 == 1))
        LRs(j_abs, level) = 1;
    else
        LRs(j_abs, level) = (L1*L2 + 1)/(L1+L2);
    end
else
    if(u(j-1) == 0)
        LRs(j_abs, level) = L2 * L1;
    else
        LRs(j_abs, level) = L2 / L1;
    end
end

L = LRs(j_abs, level);

end

% Decide according to the lr
function decoded_bit = decide(current_lr)
if(current_lr == 0)
    decoded_bit = 1;
elseif(current_lr == Inf)
    decoded_bit = 0;
elseif(current_lr == 1)
    decoded_bit = nan;
else
    error('Unexpected likelihood ratio')
end
end
