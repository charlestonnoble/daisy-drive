function [T, freqs, parms] = ...
    simulate_standard_multiple_islands(n, N, d0, P, d_cost, r_cost, ...
    t_max, r, R, q)

if nargin == 0
    
    % Drive/simulation parameters
    d0          = 0.15;
    g           = 1;
    n           = 5;        % guides
    N           = 5;        % islands
    d_cost      = 0.1;
    r_cost      = 1.0;
    q           = 1;
    P           = 0.80;
    t_max       = 500;
    plot_bool = 1;
    
    % Population size array
    R = ones(N,1);

    % Migration rate matrix
    r = 1e-2;
    r = r * (gallery('tridiag',ones(1,N-1),ones(1,N),ones(1,N-1))-eye(N));
    
else
    g = 1;
    plot_bool = 0;
end


% Set up various arrays that will be used later
a_names = gen_allele_lookup_arr(n);
d_names = gen_genotp_lookup_arr(a_names);               %#ok<*NASGU>
[parr, parr_names] = gen_p_arr(n, q, P, g, a_names);
[farr, farr_names] = gen_f_arr(n, d_cost, r_cost, a_names);
init_arr = gen_init_arr(d0,n,length(farr),N);

% Put everything into a struct that we can easily pass around
parms = struct();
parms.n = n;
parms.farr = farr;
parms.parr = parr;
parms.N = N;
parms.r = r;
parms.R = R;
parms.gt_names = d_names;
parms.ht_names = a_names;

% Run the simulation
options = odeset('RelTol', 1E-8, 'AbsTol', 1E-8);
[T,Y] = ode45(@(t, y) update_function(t, y, parms), [0 t_max], init_arr, options);

freqs = convert_to_allele_frequencies(Y, parms);

if plot_bool
    plot_sims(T,freqs,N);
end

end

function [] = plot_sims(T,freqs,islands)

addpath('Plotting Functions/')

temp = brewermap(3,'Reds');
D_clr = temp(3,:);
temp = brewermap(3,'Blues');
R_clr = temp(3,:);

figure('position',[360   449   872   249],'color','w');
ha = tight_subplot(1,islands,0.04,[0.15,0.08],[0.06,0.02]);

for ax = 1:islands
    
    pl_cell = [];
    str_cell = {};
    
    axes(ha(ax))
    hold on
    
    plot(gca,T,squeeze(freqs(1,ax,:)),'-','Color',D_clr);
    plot(gca,T,squeeze(freqs(end,ax,:)),'-','Color',R_clr);
    
    xlabel('Generations')
    set(gca,'yscale','linear')
    ylim([0,1])
    xlim([0,max(T)])
    
    if ax == 1
        ylabel('Allele frequency')
    end
    if ax == islands
        legend('Drive','Resistance')
    end
end

end

function dy = update_function(~, y, parms)

%Load parameters
n       = parms.n;
farr    = parms.farr;
parr    = parms.parr;
N       = parms.N;
R       = parms.R;
r       = parms.r;
gt_count = size(parr,1);
ht_count = size(parr,2);

y_arr = reshape(y,[gt_count,N]);

%Make the F_array, F_i stores the F_I where anames(i) is the allele name
g_arr = zeros(ht_count,N);

for j = 1:N
    y_temp = y_arr(:,j);
    g_arr(:,j) = ((y_temp.*farr)'*parr)';
end

% Density dependence
psi_arr = sum(g_arr);

% Migration terms
mig_arr = zeros(size(y_arr));
for i = 1:N
    denominator = R(i);
    for j = 1:N
        if i ~= j
            r_val = r(i,j);
            mig_arr(:,i) = mig_arr(:,i) + ...
                r_val / denominator * (y_arr(:,j)-y_arr(:,i));
        end
    end
end

dy_arr = zeros(size(y_arr));

for isl = 1:N
    idx = 1;
    for i = 1:ht_count
        for j = i:ht_count
            kd = (i == j);
            dy_arr(idx,isl) = (2-kd) * g_arr(i,isl) * g_arr(j,isl);
            idx = idx+1;
        end
    end
end
dy_arr = dy_arr + mig_arr - y_arr .* repmat(psi_arr, gt_count, 1) .^2;
dy = reshape(dy_arr,[N*gt_count,1]);

end

function init_arr = gen_init_arr(d0,n,genotype_count,N)

init_arr = zeros(genotype_count*N,1);

dr_idx1 = 1;
wt_idx1 = 2*n+3;
dr_idxs = dr_idx1:genotype_count:length(init_arr);
wt_idxs = wt_idx1:genotype_count:length(init_arr);

init_arr(wt_idxs) = 1;
init_arr(dr_idxs(1)) = d0;
init_arr(wt_idxs(1)) = 1-d0;

end

function [farr, farr_names] = gen_f_arr(n, d_cost, r_cost, alr)

f = struct();
f.aa = 1;
f.dd = 1-d_cost;
f.rr = 1-r_cost;
f.ad = 1-d_cost;
f.ar = 1;
f.dr = 1-d_cost;

allele_count = 2*n+2;
total_count = (1/2)*allele_count*(allele_count+1);

farr = zeros(total_count,1);
farr_names = cell(total_count,1);

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        str = [alr{i} alr{j}];
        farr_names{idx} = str;

        if strcmp(alr{i}(1),'D')
            if strcmp(alr{j}(1),'D')
                farr(idx) = f.dd;
            elseif strcmp(alr{j}(1),'S')
                farr(idx) = f.ad;
            elseif strcmp(alr{j}(1),'R')
                farr(idx) = f.dr;
            end
        elseif strcmp(alr{i}(1),'S')
            if strcmp(alr{j}(1),'S')
                farr(idx) = f.aa;
            elseif strcmp(alr{j}(1),'R')
                farr(idx) = f.ar;
            end
        elseif strcmp(alr{i}(1),'R')
            if strcmp(alr{j}(1),'R')
                farr(idx) = f.rr;
            end
        end
        idx = idx + 1;
    end
end

end

function [parr, parr_names] = gen_p_arr(n, q, P, g, alr)

allele_count = 2*n+2;
total_count = allele_count*((1/2)*allele_count*(allele_count+1));

parr = zeros(total_count,1);
parr_names = cell(total_count,1);

p_L = @(resist_i, lost_total, cuts, n) (n-resist_i-lost_total+1) * ...
    nchoosek(lost_total-2,cuts-2)/nchoosek(n-resist_i,cuts);
p_C = @(resist_i, cuts, q, n) nchoosek(n-resist_i,cuts)*q^cuts*...
    (1-q)^(n-resist_i-cuts);

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        for k = 1:2*n+2
            str = [alr{i} alr{j} ',' alr{k}];
            parr_names{idx} = str;
            
            %p is only interesting if a drive heterozygote
            if strcmp(alr{i},'D') %D is index 1, so only have to check i
                
                if j == i 
                    %This case represents DD individual
                    parr(idx) = (k == i);
                    
                else
                    %This case is DX where X = S_i or R_i
                    
                    if strcmp(alr{j}(1), 'R')
                        
                        if strcmp(alr{k}(1), 'R')
                            % p(DR, R)
                            
                            if j == k
                                % p(DR_i, R_i)    
                                i_val = str2double(alr{k}(2:end));
                                parr(idx) = 0.5 * (1-q)^(n-i_val);
                                
                            elseif k > j
                                % p(DR_i, R_k)
                                i_val = str2double(alr{j}(2:end));
                                j_val = str2double(alr{k}(2:end));
                                
                                if j_val == i_val+1
                                    parr(idx) = 0.5*(1-P)*p_C(i_val,1,q,n);
                                elseif j_val >= i_val+2
                                    val = 0;
                                    for c = 2:j_val-i_val
                                        val = val + ...
                                            p_L(i_val,j_val-i_val,c,n) *...
                                            p_C(i_val,c,q,n);
                                    end
                                    val = val * (1-P)/2;
                                    parr(idx) = val;
                                end
                            else
                                parr(idx) = 0;
                            end
                            
                        elseif strcmp(alr{k}(1), 'S')
                            % p(DR, S)
                            parr(idx) = 0;
                            
                        elseif strcmp(alr{k}(1), 'D')
                            % p(DR, D)
                            i_val = str2double(alr{j}(2:end));
                            parr(idx) = 0.5 + 0.5*P*(1-(1-q)^(n-i_val));
                            
                        end
                        
                    elseif strcmp(alr{j}(1),'S')
                        
                        if strcmp(alr{k}(1), 'R')
                            % p(DS, R)
                            
                            i_val = str2double(alr{j}(2:end));
                            j_val = str2double(alr{k}(2:end));
                            
                            if j_val == i_val+1
                                parr(idx) = 0.5 * p_C(i_val,1,q,n)*...
                                    (1-P)*(1-g);
                            elseif j_val >= i_val+2
                                val = 0;
                                for c = 2:j_val-i_val
                                    val = val + ...
                                        p_L(i_val,j_val-i_val,c,n) * ...
                                        p_C(i_val,c,q,n);
                                end
                                val = val * (1-P)/2;
                                parr(idx) = val;
                            else
                                parr(idx) = 0; 
                            end
                            
                        elseif strcmp(alr{k}(1), 'S')
                            % p(DS, S)
                            
                            if k == j
                                i_val = str2double(alr{k}(2:end));
                                parr(idx) = 0.5*(1-q)^(n-i_val);
                                
                            elseif k>j
                                i_val = str2double(alr{j}(2:end));
                                j_val = str2double(alr{k}(2:end));
                                if j_val == i_val + 1
                                    parr(idx) = 0.5*p_C(i_val,1,q,n)*...
                                        (1-P)*g;
                                elseif j_val >= i_val + 2
                                    parr(idx) = 0;
                                end
                            else
                                parr(idx) = 0;
                                
                            end
                                
                            
                        elseif strcmp(alr{k}(1), 'D')
                            % p(DS, D)
                            
                            i_val = str2double(alr{j}(2:end));
                            parr(idx) = 0.5 + 0.5*P* ...
                                (1-(1-q)^(n-i_val));     
                        end
                    end
                end
                               
            %Otherwise just standard inheritance    
            else
                if k == j
                    parr(idx) = parr(idx) + 1/2;
                end
                if k == i
                    parr(idx) = parr(idx) + 1/2;
                end
            end
            
            idx = idx + 1;
        end
    end
end

parr = reshape(parr,2*n+2,[])';
parr_names = reshape(parr_names,2*n+2,[])';

end

function glr = gen_genotp_lookup_arr(larr)  

allele_count = length(larr);
glr = cell((1/2)*allele_count*(allele_count+1),1);
idx = 1;
for i = 1:allele_count
    for j = i:allele_count
        str = [larr{i} larr{j}];
        glr{idx} = str;
        idx = idx+1;
    end
end

end

function alr = gen_allele_lookup_arr(n)

alr = cell(2*n+2,1);
alr{1} = 'D';
for i = 2:2+n
    alr{i} = ['S' sprintf('%d', i-2)];
end
for i = 3+n:2+2*n
    alr{i} = ['R' sprintf('%d', i-2-n)];
end

end

function out_freqs = convert_to_allele_frequencies(Y, parms)

ht_count = size(parms.parr,2);
gt_count = size(parms.parr,1);
islands = length(parms.R);

Y = permute(Y,[2,1]);
T = size(Y,2);
Y_arr = reshape(Y,[gt_count,islands,T]);
out_freqs = zeros(ht_count, islands, T);

% Make a quick lookup table
tbl = zeros(nchoosek(ht_count,2)+ht_count,2);
idx = 1;
for i = 1:ht_count
    for j = i:ht_count
        tbl(idx,1) = i;
        tbl(idx,2) = j;
        idx = idx + 1;
    end
end

for isl = 1:islands
    idx = 1;
    for i = 1:ht_count
        for j = i:ht_count
            ht1 = tbl(idx,1);
            ht2 = tbl(idx,2);
            if i == j
                out_freqs(ht1,isl,:) = out_freqs(ht1,isl,:) + Y_arr(idx,isl,:);
            else
                out_freqs(ht1,isl,:) = out_freqs(ht1,isl,:) + 0.5*Y_arr(idx,isl,:);
                out_freqs(ht2,isl,:) = out_freqs(ht2,isl,:) + 0.5*Y_arr(idx,isl,:); 
            end
            idx = idx + 1;
        end
    end
end

end