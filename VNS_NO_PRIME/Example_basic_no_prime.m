%% =====CONSTRUCT A CONCATENATED DESIGN====================================
% This script shows how to construct an concatenated designs from a regular  
% design with a number of basic factors that is a prime.
%
% WARNING!: TO BE RUN ON MATLAB 2017 OR ABOVE.
%
%% ====Example: Construct a concatenated design with 768 runs and 24 factors

% Set parameters.----------------------------------------------------------
nfactors = 24;      % Number of factors.
ngenerators = 16;   % Number of generators of the regular design.
p = 7;              % Number of factors to permute (MUST be a prime).
d = 3;              % Number of parent designs.
maxiter = 100;      % Maximum number of iterations of the VNS algorithm.
if_parallel = 0;    % Parallel computations for the algorithm: 0:No,1:Yes.
rng(442);           % Set seed.

% ====================CONSTRUCT DESIGN=====================================
% Load regular design.-----------------------------------------------------
m = num2str(nfactors);
k = num2str(ngenerators);
myfile = strcat('regular_designs/MA_m',m,'_k',k, '_d1.txt');
regular_design = textread(myfile);
b = nfactors - ngenerators;
n = 2^b;
basic_vec = 1:b;

% ====Execute VNS algorithm================================================
Bset = nchoosek(basic_vec, p);
n_bset = nchoosek(b, p);
designs = zeros(n*d, nfactors, p);
evaluate_designs = zeros(p, n*d/16 + 3);

% Test for all subsets of basic factors.-----------------------------------
tic;
for ii = 1:n_bset
    designs(:, :, ii)  = VNS_noprime( regular_design, Bset(ii,:), d, maxiter, if_parallel);
    b_remove  = basic_vec(~ismember(basic_vec, Bset(ii,:)));
    % Evaluate design.-----------------------------------------------------
    cfv = F4(designs(:, :, ii)); 
    Ffour = cfv{1}';
    evaluate_designs(ii,:) = [b_remove, Ffour(2,:), cfv{2}, rankX2(designs(:, :, ii))];
end
time = toc;

% ======================REPORT DESIGNS=====================================
disp('Size of designs');disp(size(designs));
smatFfour = sortrows(evaluate_designs, 2:size(evaluate_designs,2));
disp('F4 B4 df'); disp(smatFfour);
disp('Computing time'); disp(time);
