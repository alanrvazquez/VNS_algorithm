function [ combined ] = concatenate( designs, n, m, d )
% Concatenate several designs. 
%
% INPUTS:  
% n         Run size of the designs.
% m         Number of factors. 
% d         Number of designs.
% designs   The n-by-m-by-k array containing the designs to be concatenated.
%
% OUTPUTS:
% combined  The (n*k)-by-m concatenated design.
%==========================================================================
combined = zeros(n, m);
for ii = 1:d
    bb = n*(ii-1)+1;
    ee = n*ii;
    combined(bb:ee, :) = designs(:, :, ii);
end
end

