% Encode the input vector, by adding the frozen bits and performing the
% polar transform.
%   
% encoded_input = encode_input(input, frozen_bits, A, A_c)
% 
% INPUT
%   input           1 x K vector (in increasing order)
%   frozen_bits     1 x (BLOCKLENGTH-K) vector (in increasing order)
%   A               1 x BLOCKLENGTH logical vector
%   A_c             1 x BLOCKLENGTH logical vector 
% 
% OUTPUT
%   encoded_input   1 x BLOCKLENGTH vector
function encoded_input = encode_input(input, frozen_bits, A, A_c)
BLOCKLENGTH = length(A);

% Concatenate input and frozen bits
bits_to_combine = zeros(1, BLOCKLENGTH);
bits_to_combine(A) = input;
bits_to_combine(A_c) = frozen_bits;

% Recursively combine the bits using polar transformation
encoded_input = combine_bits(bits_to_combine, BLOCKLENGTH);

end


% Recursively combine the bits using polar transformation
% See Arikan paper figure 3
% 
% INPUT
%   u               1 x BLOCKLENGTH vector
%   BLOCKLENGTH     scalar
% 
% OUTPUT
%   x               1 x BLOCKLENGTH vector
function x = combine_bits(u, BLOCKLENGTH)
if(BLOCKLENGTH == 1)
    x = u;
    return;
end

u_odd = u(1:2:BLOCKLENGTH-1);
u_even = u(2:2:BLOCKLENGTH);

s_odd = mod(u_odd+u_even, 2);
s_even = u_even;

% Reverse shuffle operation R_N
v_first_half = s_odd;
v_second_half = s_even;

% Recursively encode v
x_first_half = combine_bits(v_first_half, BLOCKLENGTH / 2);
x_second_half = combine_bits(v_second_half, BLOCKLENGTH / 2);

x = [x_first_half x_second_half];
end
