function [result] = F4(design) 
% Compute the F4 vector and B4 value of a two-level design of strength 3.
%
% INPUTS:
% design    The N-by-m upper design with coded levels -1 and +1 (matrix).
%
% OUTPUTS:
% result    A 2-by-1 cell array. The first element of the array contains 
%           the F4 vector while the second element the B4 value. 
%
% Set up inputs and variables.---------------------------------------------
[N, m] = size(design);
nelF4=ceil(N/16);
result = cell(2, 1);
possibJ = N - 16*(0:(nelF4-1)); % Possible values for the J4-characteristics.

% Compute J4-characteristics.----------------------------------------------
B4 = 0;
F4 = zeros(1, nelF4);
nchfour = nchoosek(1:m, 4);
ncomb = nchoosek(m, 4);
J4 = zeros(1, ncomb);
for ii = 1:ncomb
    J4(ii) = abs( sum(design(:, nchfour(ii, 1)).*design(:, nchfour(ii, 2)).*design(:, nchfour(ii, 3)).*design(:, nchfour(ii, 4)) ) );
    B4 = B4 + ((sum(design(:, nchfour(ii, 1)).*design(:, nchfour(ii, 2)).*design(:, nchfour(ii, 3)).*design(:, nchfour(ii, 4)) ) )/N)^2;
end

% Compute confounding frequency vector of length 4.------------------------
for clind=1:nelF4
    F4(clind)=sum(sum(J4==N-(clind-1)*16));
end

% Svae results.------------------------------------------------------------
result{1} = [possibJ', F4'];
result{2} = B4;
end