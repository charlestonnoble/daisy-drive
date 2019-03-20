function [ ] = single_parameter_tabove_data( )

addpath('simulation functions/')

% Generate the default parameters
parms = gen_params;
thresh = 0.5;

% Pick a parameter that we're interested in
% p_string = 's';
% p_range = [-3, -1];
% log_lin = 'log';

% p_string = 'H';
% p_range = [0.7, 1];
% log_lin = 'linear';

% p_string = 'd';
% p_range = [-4, 0];
% log_lin = 'log';

% p_string = 'd0';
% p_range = [-2, 0];
% log_lin = 'log';

p_string = 'c';
p_range = [-4, 0];
log_lin = 'log';

% Find the value of the parameter where it first hits 99% frequency
[parm_thresh_val, dir] = find_thresh_val(p_string, p_range, parms, thresh, log_lin);

% Explore
[x_arr, y_arr] = find_tabove_time_multiple(parms, p_string, ...
    parm_thresh_val, p_range, thresh, dir, log_lin);

close all
plot(x_arr,y_arr,'ko-')
xlabel(p_string)
ylabel(['Generations above ' num2str(thresh)])
set(gca,'yscale','log')
set(gca,'xscale',log_lin)

save(['./time_above_50pct_figs/' p_string '.mat'])

end

function [x_arr, y_arr] = ...
    find_tabove_time_multiple(parms, parm_string, parm_thresh_val, ...
    parm_range, thresh, dir, lin_log)

pts = 400;
val1 = parm_thresh_val;
if dir == 1
    if strcmp(lin_log,'linear')
        x_arr = linspace(val1,min([100*val1, parm_range(2)]),pts);
    else
        x_arr = logspace(log10(val1),min([log10(val1)+2, parm_range(2)]),pts);
    end
else
    if strcmp(lin_log,'linear')
        x_arr = linspace(val1/1000,val1,pts);
    else
        x_arr = logspace(log10(val1)-3,log10(val1),pts);
    end
    
end
y_arr = zeros(size(x_arr));
for i = 1:length(x_arr)
    y_arr(i) = ...
        find_tabove_time_single(setfield(parms,parm_string,x_arr(i)),thresh);
    if sum(isnan(y_arr)) > 2
        y_arr(i:end) = nan;
        break
    end
    [x_arr(i),y_arr(i)]
end

end

function rtn_t = find_tabove_time_single(parms, thresh)

freq_end = inf;
len_arr = [100, 500, 1000, 5000, 10000, 50000];
for i = 1:length(len_arr)
    [T,D,~] = run_single_sim(parms, len_arr(i), 0);
    
    freqs = D(:,end);
    i1 = find(freqs > thresh, 1, 'first');
    i2 = find(freqs > thresh, 1, 'last');
    
    if ~isempty(i1) && ~isempty(i2) && D(end) < thresh && i1 ~= i2
        rtn_t = T(i2) - T(i1); 
        return
    end
    
end
rtn_t = NaN;

end

function [rtn_val, dir] = find_thresh_val(p_str, p_rge, parms, thresh, log_lin)

% First, figure out in which direction the max frequency increases
if strcmp(log_lin,'linear')
    vals = linspace(p_rge(1),p_rge(2),4);
    vals = vals(2:3);
elseif strcmp(log_lin,'log')
    vals = logspace(p_rge(1),p_rge(2),4);
    vals = vals(2:3);
end

[~,D1,~] = run_single_sim(setfield(parms,p_str,vals(1)), 1000, 1);
[~,D2,~] = run_single_sim(setfield(parms,p_str,vals(2)), 1000, 1);
max_val1 = max(D1(:,end));
max_val2 = max(D2(:,end));
if max_val2 > max_val1
    dir = 1;
elseif max_val2 < max_val1
    dir = -1;
else
    error('No change.')
end

tol = 1e-8;
if strcmp(log_lin,'log')
    val_l = 10^p_rge(1);
    val_h = 10^p_rge(2);
elseif strcmp(log_lin,'linear')
    val_l = p_rge(1);
    val_h = p_rge(2);
end
curr_diff = abs(val_h-val_l);

while curr_diff > tol
    val = (val_h+val_l)/2;
    parms = setfield(parms,p_str,val);
    [~,D,~] = run_single_sim(parms, 1000, 1);
    if max(D(:,end)) > thresh
        if dir == 1
            val_h = val;
        else
            val_l = val;
        end
    elseif max(D(:,end)) < thresh
        if dir == 1
            val_l = val;
        else
            val_h = val;
        end
    end
    curr_diff = abs(val_h - val_l);
    disp(['Current value: ' num2str((val_h+val_l)/2)])
end

rtn_val = mean([val_l val_h]);

end

function parms = gen_params()

% Uncomment to study effects of "s" cargo resistance cost
% parms.n  = 3;
% parms.d0 = 0.01;
% parms.H  = 0.95;
% parms.c  = 0.05;
% parms.d  = 0.01e-2;
% parms.s  = 0;

% Uncomment for all other parameters
parms.n  = 3;
parms.d0 = 0.1;
parms.H  = 0.95;
parms.c  = 0.05;
parms.d  = 0.01e-2;
parms.s  = 0;

end

function [T,D,R] = run_single_sim(parms, gens, cutoff)

if nargin == 1; gens = 200; cutoff = 0; end;

[T,D,R] = simulate_resistance(parms.n, parms.d0, parms.H, parms.c, parms.d, parms.s, ...
    gens, cutoff, 0);

end