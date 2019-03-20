function [T, D_freqs, R_freqs] = ...
    simulate_resistance_multiple_islands(n,N,d0,h,c,d,s,t_max,r,R)

if nargin == 0
    % Drive/simulation parameters
    n       = 3;        % daisy elements
    N       = 5;        % islands
    d0      = 0.15;      % release frequency on island 1
    h       = 0.80;      % homing efficiency
    c       = 0.1;      % cargo cost
    d       = 0.01e-2;  % upstream cost
    s       = 1;        % cargo resistance cost
    t_max   = 200;      % generations
    plot_bool = 1;      % plot yes or no
    
    % Population size array
    R = ones(N,1);

    % Migration rate matrix
    r = 1e-2;
    r = r * (gallery('tridiag',ones(1,N-1),ones(1,N),ones(1,N-1))-eye(N));
else
    plot_bool = 0;
end

% Set up various arrays that will be used later
haplotypes = init_haplotypes(n);
genotypes = init_genotypes(n);
mat = init_prod_mat(n,h,haplotypes,genotypes);
fit_arr = init_fit_arr(genotypes,c,d,s,n);

% Put everything into a struct that we can easily pass around
parms = struct();
parms.n     = n;
parms.c     = c;
parms.H     = h;
parms.R     = R;
parms.r     = r;
parms.genotypes = genotypes;
parms.mat = mat;
parms.fit_arr = fit_arr;
parms.haplotypes = haplotypes;

% Make the initial frequency array
init_arr = gen_init_arr(n, d0, R);

% Run the simulation
tol = 1E-8;
options = odeset('RelTol',tol,'AbsTol',tol);
[T,Y] = ode45(@(t, y) update_function(t, y, parms), [0 t_max], init_arr, options);

% Convert the genotype frequencies to drive allele frequencies
[D_freqs, R_freqs] = convert_to_allele_frequencies(n, genotypes, Y, R);

if plot_bool
    plot_sims(T,D_freqs,R_freqs,N);
end

end

function [] = plot_sims(T,D_freqs,R_freqs,islands)

addpath('Plotting Functions/')

n_D = size(D_freqs,1);
n_R = size(R_freqs,1);

D_colors = brewermap(n_D,'Reds');
R_colors = brewermap(n_R,'Blues');

close all; figure('position',[360   449   872   249],'color','w');
ha = tight_subplot(1,islands,0.04,[0.15,0.08],[0.06,0.02]);

for ax = 1:islands
    
    pl_cell = [];
    str_cell = {};
    
    axes(ha(ax))
    hold on
    
    for i = 1:n_D
        pl_cell(end+1) = plot(gca,T,squeeze(D_freqs(i,ax,:)),'-','Color',D_colors(i,:));
        str_cell{end+1} = ['Drive allele ' num2str(i)];
    end
    
    for i = 1:n_R
        pl_cell(end+1) = plot(gca,T,squeeze(R_freqs(i,ax,:)),'-','Color',R_colors(i,:));
        str_cell{end+1} = ['Resistant allele ' num2str(i)];
    end
    
    xlabel('Generations')
    set(gca,'yscale','linear')
    ylim([0,1])
    xlim([0,max(T)])
    
    if ax == 1
        ylabel('Allele frequency')
    end
    if ax == islands
        legend(pl_cell([n_D:-1:1, end:-1:n_D+1]), ...
            str_cell{[n_D:-1:1, end:-1:n_D+1]})
    end
end

end

function mat = init_prod_mat(n,h,haplotypes,genotypes)

mat = zeros(6^n,3^n);

for i = 1:6^n
    gt = genotypes(i,:);
    [alpha,beta] = split_genotype(gt);
    for j = 1:3^n
        b = haplotypes(j,:);
        prob = calc_prob(alpha,beta,b,h);
        mat(i,j) = prob;
    end
end

end

function gt = combine_haplotypes(hap1, hap2)
% Combines two haplotypes to produce the corresponding genotype.

gt = zeros(size(hap1));
for i = 1:length(hap1)
    gt(i) = combine_single_haplotype_pair(hap1(i),hap2(i));
end

end

function gt = combine_single_haplotype_pair(h1, h2)

h1_temp = min([h1,h2]);
h2_temp = max([h1,h2]);
h1=h1_temp;
h2=h2_temp;
if h1 == 0
    if h2 == 0
        gt = 0;
    elseif h2 == 1
        gt = 1;
    elseif h2 == 2
        gt = 2;
    end
elseif h1 == 1
    if h2 == 1
        gt = 3;
    elseif h2 == 2
        gt = 4;
    end
elseif h1 == 2
    gt = 5;
end

end

function hts = split_single_genotype(gt)

% Define W = 0; D = 1; R = 2;
switch gt
    case 0
        hts = [0, 0];
    case 1
        hts = [0, 1];
    case 2
        hts = [0, 2];
    case 3
        hts = [1, 1];
    case 4
        hts = [1, 2];
    case 5
        hts = [2, 2];
end

end

function out = calc_prob(alpha,beta,b,h)
% Calculates the genotype -> haplotype production probabilities, as
% described by Eq. (8) in the SI.

n = length(alpha);
out = 1;

% % Quickly check to see if there are any positions that would make the
% % probability zero. Specifically, if there are any WW, DD or RR positions
% % in the genotype, then only those alleles are allowed at that position.
% % We'll arrive at this conclusion via a probability of zero if we run the
% % rest of this function, but we'll save time if we just check right off the
% % bat, since most of the genotype/haplotype pairings yield p = 0.
% idxs = (alpha == beta);
% if sum(b(idxs) == alpha(idxs)) < sum(idxs)
%     out = 0;
%     return;
% end

for i = 1:n

    if i == 1
        outer_1 = 0;
        outer_2 = 1;
    else
        outer_1 = 1-g(alpha(i-1),beta(i-1),1,0);
        outer_2 =   g(alpha(i-1),beta(i-1),1,0);
    end
    
    ai = alpha(i);
    bi = beta(i);

    % Recall that 'W' = 0; 'D' = 1; 'R' = 2;
    g01 = g(ai,bi,0,1);
    g02 = g(ai,bi,0,2);
    g11 = g(ai,bi,1,1);
    g12 = g(ai,bi,1,2);
    g21 = g(ai,bi,2,1);
    g22 = g(ai,bi,2,2);
    
    first_sum = ...
                  kd(b(i),2) * g02        ...
        +         kd(b(i),2) * g21 * g01  ...
        +         kd(b(i),2) * g22        ...
        + 0.5*    kd(b(i),2) * g11 * g21  ...
        + (1-h)/2*kd(b(i),2) * g01 * g11  ...
        +         kd(b(i),1) * g12        ...
        + 0.5*    kd(b(i),1) * g11 * g21  ...
        + (1+h)/2*kd(b(i),1) * g11 * g01;
    
    second_sum = ...
              kd(b(i),0) * g02          ...
        + 0.5*kd(b(i),0) * g01 * g11    ...
        + 0.5*kd(b(i),0) * g01 * g21    ...
        +     kd(b(i),1) * g12          ...
        + 0.5*kd(b(i),1) * g11 * g01    ...
        + 0.5*kd(b(i),1) * g11 * g21    ...
        +     kd(b(i),2) * g22          ...
        + 0.5*kd(b(i),2) * g01 * g21    ...
        + 0.5*kd(b(i),2) * g11 * g21;
    
    sum_val = outer_1 * first_sum + outer_2 * second_sum;
    out = out * sum_val;
     
end

end

function out = g(a,b,c,k)
% Calculates the value of the gamma in the SI Eq. (8).

switch k
    case 0
        out = (1-(a==c))*(1-(b==c));
    case 1
        out = ((a==c)*(1-(b==c))+(b==c)*(1-(a==c)));
    case 2
        out = (a==c)*(b==c);
end

end

function out = kd(val1, val2)
% Kronecker delta syntax.

out = (val1 == val2);

end

function haplotypes = init_haplotypes(n)
% Produces an array holding all the haplotypes.

temp = 0:3^n-1;
haplotypes=dec2base(temp, 3)-'0';

end

function genotypes = init_genotypes(n)
% Produces an array holding all the genotypes.

temp = 0:6^n-1;
genotypes=dec2base(temp, 6)-'0';

end

function idx = genot_idx(hap1, hap2)

gt = combine_haplotypes(hap1, hap2);
idx = sum(gt .* 6.^((length(gt)-1):-1:0))+1;

end

function arr = init_fit_arr(genotypes,c,d,s,n)

F_vals = zeros(n,1);
F_vals(1:n-1) = 1-d;
F_vals(end) = 1-c;

K_vals = zeros(n,1);
K_vals(1:n-1) = 1;
K_vals(end) = 1-s;

arr = zeros(size(genotypes,1),1);
for g = 1:size(genotypes,1)
    gt = genotypes(g,:);
    [alpha,beta] = split_genotype(gt);
    fitness = 1;
    for i = 1:n
        ai = alpha(i);
        bi = beta(i);
        exp1 = kd(ai,1)+kd(bi,1)*(1-(ai==bi));
        exp2 = kd(ai,2)+kd(bi,2)*(1-(ai==bi));
        fitness = fitness * F_vals(i)^exp1 * K_vals(i)^exp2;
    end
    arr(g) = fitness;
end

end

function init_arr = gen_init_arr(n, drive_init, R)

init_arr = zeros(length(R)*6^n,1);
didx1 = sum(3*(6.^(0:n-1)))+1;

wt_homozygote_idxs = 1:6^n:length(init_arr);
dr_homozygote_idxs = didx1:6^n:length(init_arr);

init_arr(wt_homozygote_idxs) = 1;
init_arr(dr_homozygote_idxs(1)) = drive_init;
init_arr(wt_homozygote_idxs(1)) = 1 - drive_init;

end

function [alpha, beta] = split_genotype(genotype)
% Splits a genotype into its two constituent haplotypes.

n = size(genotype,2);
gt_arr = zeros(2,n);

for i = 1:n
    gt = genotype(i);
    hts = split_single_genotype(gt);
    gt_arr(:,i) = hts';
end

alpha = gt_arr(1,:);
beta = gt_arr(2,:);

end

function dy = update_function(~, y, parms)

% Load parameters
n               = parms.n;
R               = parms.R;
r               = parms.r;
mat             = parms.mat;
fit_arr         = parms.fit_arr;
haplotypes      = parms.haplotypes;

isl_count = length(R);
y_arr = reshape(y,[6^n,isl_count]);

% Gamete production (F^i and F^m) values
g_arr = zeros(3^n,isl_count);

for j = 1:isl_count
    y_temp = y_arr(:,j);
    g_arr(:,j) = ((y_temp.*fit_arr)'*mat)';
end

% Density dependence
psi_arr = sum(g_arr);

% Migration terms
mig_arr = zeros(size(y_arr));
for i = 1:isl_count
    denominator = R(i);
    for j = 1:isl_count
        if i ~= j
            r_val = r(i,j);
            mig_arr(:,i) = mig_arr(:,i) + ...
                r_val / denominator * (y_arr(:,j)-y_arr(:,i));
        end
    end
end

% Change in genotype frequencies
dy_arr = zeros(size(y_arr));
for isl = 1:isl_count
    for i = 1:3^n
        for j = i:3^n
            genotype_idx = genot_idx(haplotypes(i,:), haplotypes(j,:));
            if i == j
                dy_arr(genotype_idx,isl) = dy_arr(genotype_idx,isl) + g_arr(i,isl) * g_arr(j,isl);
            else
                dy_arr(genotype_idx,isl) = dy_arr(genotype_idx,isl) + 2 * g_arr(i,isl) * g_arr(j,isl);
            end
        end
    end
end

dy_arr = dy_arr + mig_arr - y_arr .* repmat(psi_arr,6^n,1) .^ 2;
dy = reshape(dy_arr,[isl_count*6^n,1]);

end

function [out_D, out_R] = convert_to_allele_frequencies(n, genotypes, Y, R)

Y = permute(Y,[2,1]);
islands = length(R);
T = size(Y,2);
Y_arr = reshape(Y,[6^n,islands,T]);

out_D = zeros(n, islands, T);
for i = 1:6^n
    for j = 1:n
        val = genotypes(i,j);
        if val == 1 || val == 4
            mult = 0.5;
        elseif val == 3
            mult = 1;
        else
            mult = 0;
        end
        if mult > 0
            for isl = 1:islands
                out_D(j,isl,:) = out_D(j,isl,:) + mult * Y_arr(i,isl,:);
            end
        end
    end
end

out_R = zeros(n, islands, T);
for i = 1:6^n
    for j = 1:n
        val = genotypes(i,j);
        if val == 2 || val == 4
            mult = 0.5;
        elseif val == 5
            mult = 1;
        else
            mult = 0;
        end
        if mult > 0
            for isl = 1:islands
                out_R(j,isl,:) = out_R(j,isl,:) + mult * Y_arr(i,isl,:);
            end
        end
    end
end

end