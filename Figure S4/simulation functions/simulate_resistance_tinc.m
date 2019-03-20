function [t90] = ...
    simulate_resistance_tinc(n,d0,h,c,d,s,t_max,cutoff_bool,cont_rel)

% Note that, throughout, we refer to alleles as 'W' = 0; 'D' = 1; 'R' = 2.
if nargin == 0
    
    % Default parameters
    n           = 4;        % number of elements
    d0          = 0.01;     % release frequency
    h           = 0.80;     % drive efficiency
    c           = 0.20;     % payload fitness cost (dominant)
    d           = 1E-2;     % upstream element fitness cost (dominant)
    s           = 0.10;     % payload-resistance cost (dominant)
    t_max       = 100;      % number of generations
    plot_bool   = 1;        % if 1, sim plotted; if 0, not plotted
    cutoff_bool = 0;        % if 1, sim stopped after any of three events:
                            %    (i) payload element begins declining
                            %   (ii) payload element exceeds 0.999
                            %  (iii) payload element drops below 1E-4
                            % if 0, simulation goes all the way to t_max
    cont_rel    = 0;        % frequency of drive to release in each gen.
                            % (if single release, set to zero)
else
    % If parameters specified, assume we don't want to plot the result here
    plot_bool = 0;
end

% If everything except continuous release specified, make it 0
% (assuming single release by default)
if nargin == 8 
    cont_rel = 0;
end

% Generate structures 
genotypes   = init_genotypes(n);
haplotypes  = init_haplotypes(n);
fit_arr     = init_fit_arr(genotypes,c,d,s,n);
mat         = init_prod_mat(n,h,haplotypes,genotypes);
[init_arr, cont_rel_arr] = gen_init_arr(n,d0,genotypes,cont_rel);
gen_order   = gen_genotype_order_arr(haplotypes);

% Create a struct to pass the parameters around
parms               = struct();
parms.genotypes     = genotypes;
parms.haplotypes    = haplotypes;
parms.gen_order     = gen_order;
parms.fit_arr       = fit_arr;
parms.mat           = mat;
parms.n             = n;
parms.cont_rel_arr  = cont_rel_arr;
parms.cont_rel_sum  = sum(cont_rel_arr);

% Run the simulation
if cutoff_bool
    options = odeset('RelTol',1E-8,'AbsTol',1E-8, 'Events', @(t,y) event(t,y,parms));
else
    options = odeset('RelTol',1E-8,'AbsTol',1E-8);
end
[T,Y] = ode45(@(t, y) update_function(t, y, parms), [0 t_max], init_arr, options);

% Convert the genotype frequencies to allele frequencies
[D_allele_freqs,R_allele_freqs] = convert_to_allele_frequencies(n, genotypes, Y);

[~,idx] = max(D_allele_freqs(:,end));
t90 = T(idx);

if isempty(t90)
    t90 = NaN; 
end

if plot_bool
    plot_sims(T,D_allele_freqs,R_allele_freqs(:,2:end));
end

end

function [value,isterminal,direction] = event(t,Y,parms)

genotypes = parms.genotypes;
n = parms.n;
dY = update_function(t,Y,parms);

idxsHet = genotypes(:,n)==1 | genotypes(:,n)==4;
idxsHom = genotypes(:,n)==3;
NDriveFreqChange = (sum(dY(idxsHet)) + 2 * sum(dY(idxsHom))) / 2;
NDriveFreq = (sum(Y(idxsHet)) + 2 * sum(Y(idxsHom))) / 2;

value = [NDriveFreqChange, (1-1E-3)-NDriveFreq, (1E-4)-NDriveFreq];
isterminal = [1, 1, 1];
direction = [-1, -1, 1];

end

function [init_arr,cont_rel_arr] = gen_init_arr(n,d0,genotypes,cont_rel)

init_arr = zeros(6^n,1);
init_arr(1) = 1-d0;
all_drive_idx = sum(3*ones(6^n,n)==genotypes,2) == n;
init_arr(all_drive_idx) = d0;
cont_rel_arr = zeros(6^n,1);
cont_rel_arr(all_drive_idx) = cont_rel;

end

function [] = plot_sims(T,D,R)

n_D = size(D,2);
n_R = size(R,2);

close all;
figure;
hold on;
ax = gca;

pl_cell = [];
str_cell = {};
for i = 1:n_D
    pl_cell(end+1) = plot(ax,T,D(:,i),'-');
    str_cell{end+1} = ['Drive allele ' num2str(i)];
end
for i = 1:n_R
    pl_cell(end+1) = plot(ax,T,R(:,i),'-');
    str_cell{end+1} = ['Resistant allele ' num2str(i+1)];
end
ylabel('Allele frequency')
xlabel('Generations')
leg=legend(pl_cell([n_D:-1:1, end:-1:n_D+1]),str_cell{[n_D:-1:1, end:-1:n_D+1]});
set(leg,'location','southwest')
set(gca,'yscale','log')

end

function dy = update_function(~, y, parms)

% Load parameters
n               = parms.n;
fit_arr         = parms.fit_arr;
mat             = parms.mat;
gen_order       = parms.gen_order;
cont_rel_arr    = parms.cont_rel_arr;
cont_rel_sum    = parms.cont_rel_sum;

% Gamete production (g) values
g_arr = zeros(1,3^n);
for i = 1:6^n       % Go over all genotypes
    g_arr = g_arr + y(i) * fit_arr(i) * mat(i,:);
end
psi = sum(g_arr);

% Change in genotype frequencies
dy = zeros(size(y));
idx = 1;
for i = 1:3^n
    for j = i:3^n
        genotype_idx = gen_order(idx);
        if i == j
            dy(genotype_idx) = dy(genotype_idx) + g_arr(i) * g_arr(j);
        else
            dy(genotype_idx) = dy(genotype_idx) + 2 * g_arr(i) * g_arr(j);
        end
        idx = idx+1;
    end
end
dy = dy + cont_rel_arr;
dy = dy - (psi^2 + cont_rel_sum) * y;

end

function order_arr = gen_genotype_order_arr(haplotypes)

n = size(haplotypes,2);
order_arr = zeros(2*3^n);
idx = 1;
for i = 1:3^n
    for j = i:3^n
        genotype_idx = genot_idx(haplotypes(i,:), haplotypes(j,:));
        order_arr(idx) = genotype_idx;
        idx = idx + 1;
    end
end

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

function genotypes = init_genotypes(n)
% Produces an array holding all the genotypes.

temp = 0:6^n-1;
genotypes=dec2base(temp, 6)-'0';

end

function haplotypes = init_haplotypes(n)
% Produces an array holding all the haplotypes.

temp = 0:3^n-1;
haplotypes=dec2base(temp, 3)-'0';

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

function [out_D, out_R] = convert_to_allele_frequencies(n, genotypes, Y)

out_D = zeros(size(Y,1), n);
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
        if val > 0
            out_D(:,j) = out_D(:,j) + mult * Y(:,i); 
        end
    end
end

out_R = zeros(size(Y,1), n);
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
        if val > 0
            out_R(:,j) = out_R(:,j) + mult * Y(:,i); 
        end
    end
end

end