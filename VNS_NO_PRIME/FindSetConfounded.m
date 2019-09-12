function [ sset ] = FindSetConfounded( design, N, nchoose4, ncomb4 )
% Find completely aliased words of length 4.  
%
% INPUTS:  
% design    The N-by-m two-level design.
% N         Number of runs. 
% ncomb4    The number of combinations of m in 4.
% nchoose4  The ncomb4-by-4 matrix of 4-factor subsets among m factors.
%
% OUTPUTS:
% sset      The 4-factor sets with a J4-characteristic equal to N.
%==========================================================================
Jk = zeros(1, ncomb4);
for ii = 1:ncomb4
    Jk_set = design(:, nchoose4(ii,1)).*design(:, nchoose4(ii,2)).*design(:, nchoose4(ii,3)).*design(:, nchoose4(ii,4)); 
    Jk(ii) = abs(sum(Jk_set));
end
sset = nchoose4(Jk == N,:);
end