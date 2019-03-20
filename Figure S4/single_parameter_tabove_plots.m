function [ ] = single_parameter_tabove_plots( )

addpath('../../Plotting Functions/export_fig/')

params = {'c','d0','d','H','s'};
param_names = {'Cargo cost','Release frequency','Upstream element cost', ...
    'Homing efficiency','Cargo resistance cost'};

[x_cell, y_cell] = load_data(params);
xlog = [1, 1, 1, 0, 1];
ylog = [1, 1, 1, 1, 1];

close all;
figure('position',[390   299   792   187],'color','w')
ha = tight_subplot(1,length(params),[0.1,0.04],[0.20,0.05],[0.05, 0.02]);

clrs = brewermap(3,'BrBG');
for i = 1:length(params)
    x = x_cell{i}; y = y_cell{i};
    axes(ha(i))
    plot(x,y,'o-','color',clrs(3,:))
    if i == 1
        ylabel('Generations above 0.5') 
    end
    xlabel(param_names{i})
    if xlog(i); set(gca,'xscale','log'); end;
    if ylog(i); set(gca,'yscale','log'); end;
    
    % Curve fits
    if i == 1
        temp = fit_1(x,y);
        hold on
        plot(x,temp(x),'k--','linewidth',1)
        leg = legend('Simulations','~1/cost');
        set(leg,'box','off')
        set(leg,'location','southwest')
    end

    set(gca,'tickdir','out')
    
end

end

function out = fit_1(x,y)

fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1]);
ft = fittype('b/x','options',fo);
out = fit(x',y',ft);

end

function [x_cell, y_cell] = load_data(params)

x_cell = cell(length(params),1);
y_cell = cell(length(params),1);

for i = 1:length(params)
    param = params{i};
    f = load(['./time_above_50pct_figs/' param '.mat']);
    x = f.x_arr;
    y = f.y_arr;
    idxs = ~isnan(y);
    x = x(idxs);
    y = y(idxs);
    x_cell{i} = x;
    y_cell{i} = y;
end

end

function parms = gen_params()

% Uncomment to study effects of "s" cargo resistance cost
parms.n  = 3;
parms.d0 = 0.01;
parms.H  = 0.95;
parms.c  = 0.05;
parms.d  = 0.01e-2;
parms.s  = 0;

% % Uncomment for all other parameters
% parms.n  = 3;
% parms.d0 = 0.1;
% parms.H  = 0.95;
% parms.c  = 0.05;
% parms.d  = 0.01e-2;
% parms.s  = 0;

end