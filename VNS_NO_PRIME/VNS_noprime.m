function [ overallbest ] = VNS_noprime( regular_design, B, d, maxiter, par_pool )
% Main function for the VNS Algorithm for improving concatenated designs.
% Version for regular designs with a number of basic factors NOT a prime.
%
% INPUTS:
% regular_design    The n-by-m two-level regular design.
% B                 Set of permuted basic factors.
% d                 Number of parent designs.
% maxiter           Number of iterations. 
% par_pool          Should the iterations of the algorithm be run in parallel?
%                   0: No, 1: Yes. 
%
% OUTPUTS:
% overallbest   The (n*d)-by-m improved design. 
%==========================================================================

%% ================== SET PARAMETERS ======================================
[n, nFactors] = size(regular_design);
N = d*n;
b = log2(size(regular_design,1));
m = nFactors - b;
gen_des = regular_design(:,(b+1):end);
nchoose4 = nchoosek(1:nFactors,4);
ncomb4 = nchoosek(nFactors,4);
% Parent designs containing only generated factors.------------------------
gen_parents = zeros(n,m,d);
for ii = 1:d
    gen_parents(:,:,ii) = gen_des;
end
cdes_generated = concatenate(gen_parents, n, size(gen_des,2), d);

% Generate basic factor design with fixed and permuted factors.------------
basic_des = regular_design(:,1:b);
basic_des_perm = basic_des(:, B);
fixed_basic = ~ismember(1:b, B);
BB = zeros(N, b);
for ii = 1:d
    bb = n*(ii-1)+1;
    ee = n*ii;
    BB(bb:ee, :) = [basic_des_perm(:, circshift( 1:length(B), ii-1)), basic_des(:,fixed_basic)];
end

% Calculate possible values for the J4-characteristics.--------------------
nvec = n*ones(1, d+1);
possibJhc = (0:d).*nvec + 1;
expvec = 5*(0:(d-1));
vecbigm = 10.^expvec; 

%% ============= CALCULATE SET OF FULLY ALIASED WORDS =====================
combined_design = [BB, cdes_generated];
sset = FindSetConfounded( combined_design, N, nchoose4, ncomb4 );
nset = size(sset,1);

%% ====================== EXECUTE ALGORITHM ===============================
concatenated_designs = zeros(N, nFactors, maxiter);
obj_vec = zeros(1, maxiter);
labsols = 1:maxiter;
objvalue = 10^40;
if par_pool % Parallel implementation.-------------------------------------
    parfor ii = 1:maxiter
        [cdesign, obj] = vnsalgorithmfour( gen_parents, BB, objvalue, n, d, m, sset, nset, possibJhc, vecbigm );
        concatenated_designs(:,:,ii) = cdesign;
        obj_vec(ii) = obj;
    end
else % Standard implementation.--------------------------------------------
    for ii = 1:maxiter 
        [cdesign, obj] = vnsalgorithmfour( gen_parents, BB, objvalue, n, d, m, sset, nset, possibJhc, vecbigm );
        concatenated_designs(:,:,ii) = cdesign;
        obj_vec(ii) = obj;
    end
end
% Save the best combined design.-------------------------------------------
sel_des = labsols(min(obj_vec) == obj_vec);
overallbest = concatenated_designs(:, :, min(sel_des));
end
