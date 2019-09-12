%% =====CONSTRUCT A CONCATENATED DESIGN====================================
% This script shows how to construct an concatenated designs from a regular  
% design with a number of basic factors that is a prime.
%
% WARNING!: TO BE RUN ON MATLAB 2017 OR ABOVE.
%
%% ====Example: Construct a concatenated design with 384 runs and 20 factors

% Set parameters.----------------------------------------------------------
nfactors = 20;      % Number of factors.
ngenerators = 13;   % Number of generators of the regular design.
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

% Isolate generated factor columns.----------------------------------------
D = regular_design(:,(b+1):end);
n = 2^b;
parents = zeros(n,ngenerators,d);
for ii = 1:d
    parents(:,:,ii) = D;
end

% Execute VNS algorithm.---------------------------------------------------
tic;
best_cdes = VNS( parents, maxiter, if_parallel);
time = toc;

% Build design.------------------------------------------------------------
basic_des = regular_design(:,1:b);
basic_parents = zeros(n, b, d);
for ii = 1:d
    basic_parents(:,:,ii) = basic_des(:, circshift( (1:b),ii-1));
end
BB = concatenate(basic_parents, n, b, d);
design = [BB, best_cdes];

% ======================REPORT DESIGN======================================
disp('Size'); disp(size(design));
evl_des = F4(design); 
disp('F4'); disp(evl_des{1});
disp('B4'); disp(evl_des{2});
disp('Rank 2FI matrix'); disp(rankX2(design));
disp('Computing time'); disp(time);
