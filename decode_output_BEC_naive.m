% Decode the received vector in a naive way (without dynamic programming).
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
function decoded_output = decode_output_BEC_naive(received_output, frozen_bits, A, A_c)

BLOCKLENGTH = length(received_output);

% Put frozen bits in a 1 x BLOCKLENGTH vector, at positions A_c, for easier
% access
frozen_bits_expanded = nan(1,BLOCKLENGTH);
frozen_bits_expanded(A_c) = frozen_bits;

% Decode bit by bit in a naive way (not reusing previous results)
decoded_output = nan(1, BLOCKLENGTH);
for j = 1:BLOCKLENGTH
    if(A_c(j))
        % If the bit is frozen, we don't need to compute anything
        decoded_output(j) = frozen_bits_expanded(j);
    else
        % To decode, first compute the likelihood ratio using the previously
        % decoded bits
        current_lr = compute_lr(received_output, decoded_output(1:j-1), BLOCKLENGTH, j);
        
        % Then decide according to the lr
        decoded_output(j) = decide(current_lr);
        
        % If we cannot recover (erasure), we stop decoding
        if(isnan(decoded_output(j)))
            break;
        end
    end
end

% Only return the information bits
decoded_output = decoded_output(A);

end


% Compute the likelihood ration using channels outputs and decoded bits
% See Arikan formula 74 and 75
%
% INPUT
%   y           1 x N vector (reduced by 2 at each recursive call)
%   u           1 x (j-1) vector
%   N           scalar
%   j           scalar
%
% OUTPUT
%   L           scalar
function L = compute_lr(y, u, N, j)

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
    
    return;
end

% Use formula 74 or 75 according to the parity of j
if(mod(j,2) == 1)
    u_odd = u(1:2:j-2);
    u_even = u(2:2:j-1);
    
    L1 = compute_lr(y(1:N/2), mod(u_odd + u_even, 2), N/2, (j+1)/2);
    
    L2 = compute_lr(y(N/2+1:N), u_even, N/2, (j+1)/2);
    
    % Use table decision for border cases (make diagram to understand)
    if((L1 == 0 && L2 == 0) || (isinf(L1) && isinf(L2)))
        L = Inf;
    elseif((L1 == 0 && isinf(L2)) || (isinf(L1) && L2 == 0))
        L = 0;
    elseif((L1 == 1 && isinf(L2)) || (isinf(L1) && L2 == 1))
        L = 1;
    else
        L = (L1*L2 + 1)/(L1+L2);
    end
else
    u_odd = u(1:2:j-3);
    u_even = u(2:2:j-2);
    
    L1 = compute_lr(y(1:N/2), mod(u_odd + u_even, 2), N/2, j/2);
    
    L2 = compute_lr(y(N/2+1:N), u_even, N/2, j/2);
    
    if(u(j-1) == 0)
        L = L2 * L1;
    else
        L = L2 / L1;
    end
end

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
