function [ ] = gen_maintext_figure( )

addpath('../../')
addpath('../../Plotting Functions/')
addpath('../../Plotting Functions/export_fig/')

gen_data = 0;

if gen_data

    % Parameters that are the same for all simulations
    N = 5;              % islands
    P = 0.80;           % homing efficiency
    c = 0.1;            % drive cost
    t_max = 500;        % simulation length
    s = 1;              % resistance cost
    R = ones(N,1);      % island sizes
    mig_rate = 10^-2;   % migration rate
    
    r = mig_rate * (gallery('tridiag',ones(1,N-1),ones(1,N),ones(1,N-1))-eye(N));
    d0 = 0.15;

    % Daisy drive
    daisy_elems = 3;
    [T_daisy, D_daisy, R_daisy] = ...
        simulate_resistance_multiple_islands(daisy_elems, N, d0, P, c, ...
        0.01e-2, s, t_max, r, R);

    % Standard drive
    [T_standard, D_standard] = simulate_standard_multiple_islands(5, N, d0, ...
        P, c, s, t_max, r, R, 1);

    % High-release no-drive
    d0_high = 0.999;
    [T_high, D_high] = simulate_standard_multiple_islands(5, N, d0_high, P, ...
        c, s, t_max, r, R, 0);
    
    save(['data_' num2str(mig_rate) '.mat'],'T_daisy','D_daisy','R_daisy','T_standard',...
        'D_standard','T_low','D_low','T_high','D_high','N','P','c',...
        't_max','s','R','r','d0','daisy_elems','d0_high')
    
else
    
    f = load('data_0.01.mat'); %#ok<*UNRCH>
    T_daisy = f.T_daisy;
    D_daisy = f.D_daisy;
    R_daisy = f.R_daisy;
    T_standard = f.T_standard;
    D_standard = f.D_standard;
    T_low = f.T_low;
    D_low = f.D_low;
    T_high = f.T_high;
    D_high = f.D_high;
    N = f.N;
    
end

close all; figure('position',[748   148   424   458],'color','w');
ha = tight_subplot(N,3,0.035,[0.08,0.02],[0.08,0.02]);

ys = 'linear';
yl = [0, 1];
lw = 1;

plot_daisy(ha(1:3:end), T_daisy, D_daisy, R_daisy, ys, yl, lw)
plot_standard(ha(2:3:end), T_standard, D_standard, ys, yl, lw)
plot_standard(ha(3:3:end), T_high, D_high, ys, yl, lw)

end

function [] = plot_standard(ha, T, freqs, ys, yl, lw)

temp = brewermap(5,'BrBG');
D_clr = temp(5,:);
R_clr = temp(1,:);

for ax = 1:size(freqs,2)
    
    pl_cell = [];
    str_cell = {};
    
    axes(ha(ax))
    hold on
    
    plot(gca,T,squeeze(freqs(1,ax,:)),'-','Color',D_clr,'linewidth',lw);
    
%     if sum(freqs(end,ax,:)) > 0
%         plot(gca,T,squeeze(freqs(end,ax,:)),'-','Color',R_clr,'linewidth',lw);
%     end
    
    set(gca,'yscale',ys)
    ylim(yl)
    xlim([0,max(T)])
    set(gca,'xtick',0:100:500)
    
    grid on
    set(gca,'yminorgrid','off')
    
    set(gca,'yticklabel',[])
    if ax == size(freqs,2)
        xlabel('Generations')
        leg=legend('D','R');
        set(leg,'box','off')
    else
        set(gca,'xticklabel',[])
    end
    
    
end

end

function [] = plot_daisy(ha, T, D_freqs, R_freqs, ys, yl,lw)

n_D = size(D_freqs,1);
n_R = size(R_freqs,1);

colors = brewermap(n_D+n_R,'BrBG');
R_colors = colors(n_D:-1:1,:);
D_colors = colors(end-n_R+1:end,:);

for ax = 1:size(D_freqs,2)
    
    pl_cell = [];
    str_cell = {};
    
    axes(ha(ax))
    hold on
    
    for i = 1:n_D
        plot(gca,T,squeeze(D_freqs(i,ax,:)),'-','Color',D_colors(i,:),'linewidth',lw);
    end
%     
%     for i = 1:n_R
%         if sum(R_freqs(i,ax,:)) > 0
%             plot(gca,T,squeeze(R_freqs(i,ax,:)),'-','Color',R_colors(i,:),'linewidth',lw);
%         end
%     end
    
    set(gca,'yscale',ys)
    ylim(yl)
    xlim([0,max(T)])
    ylabel('Frequency')
    set(gca,'xtick',0:100:500)
    
    grid on
    set(gca,'yminorgrid','off')
    
    if ax == size(D_freqs,2)
        xlabel('Generations')
        leg=legend('C','B','A','R (B)','R (A)');
        set(leg,'box','off')
    else
        set(gca,'xticklabel',[])
    end
    
end

end

