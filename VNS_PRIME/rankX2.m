function [ rango ] = rankX2( D )
% Calculate rank of the two-factor interaction matrix of a design0
%
% INPUTS:  
% D         The N-by-m two-level design.
%
% OUTPUTS:
% rango     The rank of the matrix.
%==========================================================================

m = size(D, 2); 
mod.mat = x2fx(D, 'i'); 
TWOFI = mod.mat(:, (m+2):end); 
rango = rank(TWOFI); 
end

