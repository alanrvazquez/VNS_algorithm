function [ overallbest ] = VNS( parents, maxiter, par_pool )
% Main function for the VNS Algorithm for improving concatenated designs.
% Version for regular designs with a prime number of basic factors.
%
% INPUTS:
% parents       The n-by-m-by-d array with the two-level parent designs.
% maxiter       Number of iterations. 
% par_pool      Should the iterations of the algorithm be run in parallel?
%               0: No, 1: Yes. 
%
% OUTPUTS:
% overallbest   The (n*d)-by-m improved design. 
%==========================================================================
%% ============== SET PARAMETERS ==========================================
[n,m,d] = size(parents); 
nvec = n*ones(1, d+1);
% Calculate possible values for the J4-characteristics.--------------------
possibJhc = (0:d).*nvec + 1;
expvec = 5*(0:(d-1));
vecbigm = 10.^expvec;
nchoose4 = nchoosek(1:m,4);
ncomb4 = nchoosek(m,4);
N = d*n;

%% ============= CALCULATE SET OF FULLY ALIASED WORDS =====================
cdes = concatenate(parents, n, m, d);
sset = FindSetConfounded( cdes, N, nchoose4, ncomb4 ); 
nset = size(sset,1);

%% ====================== EXECUTE ALGORITHM ===============================
concatenated_designs = zeros(N, m, maxiter);
obj_vec = zeros(1, maxiter);
labsols = 1:maxiter;
objvalue = 10^40;
if par_pool % Parallel implementation.-------------------------------------
    parfor ii = 1:maxiter
        [cdesign, obj] = vnsalgorithmfour( parents, objvalue, n, d, m, sset, nset, possibJhc, vecbigm );
        concatenated_designs(:,:,ii) = cdesign;
        obj_vec(ii) = obj;
    end
else % Standard implementation.--------------------------------------------
    for ii = 1:maxiter
        [cdesign, obj] = vnsalgorithmfour( parents, objvalue, n, d, m, sset, nset, possibJhc, vecbigm );
        concatenated_designs(:,:,ii) = cdesign;
        obj_vec(ii) = obj;
    end
end
% Save the best combined design.-------------------------------------------
sel_des = labsols(min(obj_vec) == obj_vec);
overallbest = concatenated_designs(:, :, min(sel_des));
end
