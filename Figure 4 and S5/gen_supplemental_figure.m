function [ ] = gen_supplemental_figure( )

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
    mig_rate_arr = logspace(-4,-1,13);   % migration rate
    d0 = 0.15;
    d0_high = 0.999;
    
    results_daisy = cell(length(mig_rate_arr),1);
    results_standard = cell(length(mig_rate_arr),1);
    results_inundative = cell(length(mig_rate_arr),1);
    
    for mig_index = 1:length(mig_rate_arr)
        
        mig_rate = mig_rate_arr(mig_index);
        r = mig_rate * (gallery('tridiag',ones(1,N-1),ones(1,N),ones(1,N-1))-eye(N));
        
         % Daisy drive
        daisy_elems = 3;
        [T_daisy, D_daisy, R_daisy] = ...
            simulate_resistance_multiple_islands(daisy_elems, N, d0, P, c, ...
            0.01e-2, s, t_max, r, R);
        
        results = struct();
        results.T = T_daisy;
        results.D = D_daisy;
        results.R = R_daisy;
        results_daisy{mig_index} = results;

        % Standard drive
        [T_standard, D_standard] = simulate_standard_multiple_islands(5, N, d0, ...
            P, c, s, t_max, r, R, 1);
        
        results = struct();
        results.T = T_standard;
        results.D = D_standard;
        results_standard{mig_index} = results;

        % High-release no-drive
        [T_high, D_high] = simulate_standard_multiple_islands(5, N, d0_high, P, ...
            c, s, t_max, r, R, 0);

        results = struct();
        results.T = T_high;
        results.D = D_high;
        results_inundative{mig_index} = results;
        
        % Save the results as we go
        save('data_SI.mat','results_daisy','results_standard', ...
            'results_inundative','N','P','c','t_max','s','R', ...
            'mig_rate_arr','d0','d0_high')
        
    end
    
else
    
    f = load('data_SI.mat'); %#ok<*UNRCH>
    results_daisy = f.results_daisy;
    results_standard = f.results_standard;
    results_inundative = f.results_inundative;
    mig_rate_arr = f.mig_rate_arr;
    N = f.N;
    t_max = f.t_max;
    
end

close all; figure('position',[748   148   424   458],'color','w');
ha = tight_subplot(N,3,0.035,[0.08,0.02],[0.09,0.02]);

ys = 'log';
yl = [10^-6, 1];
lw = 1;

plot_standard(ha(2:3:end), results_standard, mig_rate_arr, N, ys, yl, lw)
plot_standard(ha(3:3:end), results_inundative, mig_rate_arr, N, ys, yl, lw)
plot_daisy(ha(1:3:end), results_daisy, mig_rate_arr, N, ys, yl, lw)


end

function [] = plot_standard(ha, results, mig_arr, N, ys, yl, lw)


temp = brewermap(5,'BrBG');
clr = temp(5,:);


for ax = 1:N
    
    % Get the max-frequency curve as a function of the migration rate
    max_freq_arr = zeros(length(mig_arr),1);
    for m = 1:length(mig_arr)
        max_freq_arr(m) = max(squeeze(results{m}.D(1,ax,:)));
    end
    
    plot(ha(ax),mig_arr,max_freq_arr,'o-','Color',clr,'linewidth',lw);
    
    set(ha(ax),'xscale','log')
    set(ha(ax),'yscale',ys)
    set(ha(ax),'ylim',yl)
    set(ha(ax),'xlim',[min(mig_arr),max(mig_arr)])
    
    set(ha(ax),'xgrid','on')
    set(ha(ax),'ygrid','on')
    set(ha(ax),'yminorgrid','off')
    set(ha(ax),'xminorgrid','off')
    set(ha(ax),'yticklabel',[])
    
    if ax == N
        axes(ha(ax));
        xlabel('Migration rate')
    else
        set(ha(ax),'xticklabel',[])
    end
    
end

end

function [] = plot_daisy(ha, results, mig_arr, N, ys, yl, lw)

n_D = 3;
n_R = 3;

colors = brewermap(n_D+n_R,'BrBG');
R_colors = colors(1:n_D,:);
D_colors = colors(end-n_R+1:end,:);

for ax = 1:N
    
    axes(ha(ax)); hold on;
    
    % Get the maximum-frequency arrays
    max_freq_arr_D = zeros(length(mig_arr),3);
    max_freq_arr_R = zeros(length(mig_arr),3);
    for elem = 1:3
        for m = 1:length(mig_arr)
            max_freq_arr_D(m,elem) = max(squeeze(results{m}.D(elem,ax,:)));
            if elem > 1
                max_freq_arr_R(m,elem-1) = max(squeeze(results{m}.R(elem,ax,:)));
            end 
        end
    end

    for i = 1:n_D
        plot(gca,mig_arr,max_freq_arr_D(:,i),'o-','Color',D_colors(i,:),'linewidth',lw);
    end
    
    for i = 1:n_R-1
        plot(gca,mig_arr,max_freq_arr_R(:,i),'o--','Color',R_colors(i+1,:),'linewidth',lw);
    end
    
    set(ha(ax),'xscale','log')
    set(ha(ax),'yscale',ys)
    set(ha(ax),'ylim',yl)
    set(ha(ax),'xlim',[min(mig_arr),max(mig_arr)])
    
    set(ha(ax),'xgrid','on')
    set(ha(ax),'ygrid','on')
    set(ha(ax),'yminorgrid','off')
    set(ha(ax),'xminorgrid','off')
    
    ylabel('Max. frequency')
    
    if ax == N
        xlabel('Migration rate')
    else
        set(gca,'xticklabel',[])
    end
    
end

end

