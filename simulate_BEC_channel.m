% Simulate a BEC channel by erasing (=> NaN) each component with
% probability EPSILON.
% 
% received_output = simulate_BEC_channel(encoded_input, EPSILON)
% 
% INPUT
%   encoded_input   1 x BLOCKLENGTH vector
%   EPSILON         scalar
% 
% OUTPUT
%   received_output 1 x BLOCKLENGTH vector
function received_output = simulate_BEC_channel(encoded_input, EPSILON)
received_output = encoded_input;

received_output(rand(1, length(encoded_input)) < EPSILON) = nan;
end
