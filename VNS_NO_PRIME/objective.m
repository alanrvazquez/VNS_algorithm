function [ objval  ] = objective( design, sset, nset, possibJhc, vecbigm )
% Calculate F4 vector over selected 4-factor sets.  
%
% INPUTS:  
% design    The N-by-m two-level design.
% nset      Number of 4-factor sets. 
% sset      The nset-by-4 matrix of 4-factor sets.
% possibJhc Possible values of the J4-characteristics.
% vecbigm   The 1-by-nelF4 vector with the weight for each component of the
%           F4 vector, where nelF4 is the number of elements in the 
%           F4 vector of a strength-3 design.
%
% OUTPUTS:
% objval    The objective value.
%==========================================================================
Jk = zeros(1, nset);
for ii = 1:nset
    Jk_set = design(:, sset(ii,1)).*design(:, sset(ii,2)).*design(:, sset(ii,3)).*design(:, sset(ii,4)); 
    Jk(ii) = abs(sum(Jk_set));
end
Fkvec = histcounts(Jk, possibJhc);
objval = (vecbigm)*Fkvec';
end

