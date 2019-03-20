function [ ] = gen_fig_foldchange()

addpath('../../Plotting Functions/')
close all

% Parameters
n_arr       = [1, 2, 3, 4, 5];
final_cost  = 0.10;
upstr_cost  = 0.01e-2;
P       	= 0.60;
init_arr    = logspace(-4,0,40);
C           = init_arr;
t_max       = 20;

cm_hex = {'#c7eae5','#80cdc1','#35978f','#01665e'};
cm_rgb = cellfun(@hex2rgb,cm_hex,'uniformoutput',0);

max_arr = zeros(length(n_arr),length(init_arr));
for i = 1:length(n_arr)
    for j = 1:length(init_arr)
        [~,Y] = simulate_resistance(n_arr(i), init_arr(j), P, final_cost, upstr_cost, 0, t_max, 0, 0);
        max_arr(i,j) = Y(end,end);
        if ~mod(j,10); j, end;
    end
    i
end

figure('color','w');
ha = tight_subplot(1,1,[0,0],[0.15,0.05],[0.12,0.05]);

%
plot_vals(ha(1), init_arr, max_arr, 1, cm_rgb, '-');
%

max_arr = zeros(length(n_arr),length(init_arr));
for i = 1:length(n_arr)
    for j = 1:length(init_arr)
        [~,Y] = simulate_resistance(n_arr(i), init_arr(j), P, final_cost, upstr_cost, 0, t_max, 0, init_arr(j));
        max_arr(i,j) = Y(end,end);
        if ~mod(j,10); j, end;
    end
    i
end

%
plot_vals(ha(1), init_arr, max_arr, 0, cm_rgb, '--');
%

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
set(h,'position',[11.1528    3.9861    4.3333    3.2222])

end

function [] = plot_vals(h, init_arr, max_arr, ylab, cm, ls)

axes(h); hold on;
basevals = max_arr(1,:);
plot_arr = zeros(size(max_arr,1)-1,1);
for i = 2:size(max_arr,1)
    vals = max_arr(i,:);
    plot_arr(i-1) = plot(init_arr,vals ./ basevals, ...
        'linewidth', 1, 'color', cm{i-1}, 'linestyle', ls);
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Release frequency')
ylim([1, 3000])

if ylab
    ylabel('Fold change over single element release')
end

leg_text = cell(size(max_arr,1)-1,1);
for i=2:size(max_arr,1)
    val = [num2str(i) ' elements'];
    leg_text{i-1} = val;
end
leg = legend(plot_arr(end:-1:1), leg_text(end:-1:1));
set(leg,'box','off')

end