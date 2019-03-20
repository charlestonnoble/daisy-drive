function [] = gen_fig_2()

addpath('../../Plotting Functions/')
addpath('../..')

% Either generate or load data.
gen_data = 0;

% Parameters
n           = 3;
init_arr    = [0.02, 0.02, 0.15];
H       	= [0.95, 0.60, 0.60];
cost_d_up   = 0.01e-2; %.01 percent
cost_d_pay  = 0.08;
cost_r_end  = 1;
t_max       = 100;
phase_init  = [0.01, 0.05, 0.1];
plot_opts = {'TickDir','Out','TickLength',[0.0200 0.0250]};
leg_opts = {'FontSize',8,'Box','Off'};

D_cell = cell(3,1);
T_cell = cell(3,1);
R_cell = cell(3,1);
for i = 1:3
    [T,Y_D,Y_R] = simulate_resistance(n, init_arr(i), H(i), cost_d_pay, ...
        cost_d_up, cost_r_end, t_max, 0);
    T_cell{i} = T;
    D_cell{i} = Y_D;
    R_cell{i} = Y_R;
end

close all; figure('position',[52   243   611   357],'color','w');
ha = tight_subplot(2,3,[0.1, 0.045],[0.1,0.05],[0.07,0.1]);
cm = brewermap(4,'BrBG');
cm = cm([1,2,4],:);

for i = 1:3
    axes(ha(i)); hold on;
    T = T_cell{i};
    Y_D = D_cell{i};
    Y_R = R_cell{i};
    pl_arr = [];
    for j = 1:3
        pl_arr(end+1)=plot(T,Y_D(:,j),'color',cm(j,:),'linewidth',1);
        if j>1
            pl_arr(end+1)=plot(T,Y_R(:,j),'color',cm(j,:),'linewidth',1,'linestyle','--');
        end
        if j == 3
            disp(['Max A element frequency: '  num2str(max(Y_D(:,end)))])
        end
    end
    set(gca,'ylim',[1E-4, 1])
    set(gca,'xlim',[0,100])
    set(gca,'yscale','log')
    
    if i == 1
        ylabel('Allele frequency') 
    end
    xlabel('Time (generations)')
    
    if i == 1
        leg = legend(pl_arr,{'C','B','B resist','A','A resist'}); 
        set(leg,leg_opts{:}); 
    end
    set(gca,plot_opts{:})
end

beginH = 0.5;
steps = 200;
iterH = (1-beginH)/steps;

endval = 0.2;
iterCost = endval / steps;

if gen_data %#ok<*UNRCH>
    img1 = gen_max_freq_arr(n, 0:iterCost:endval, beginH:iterH:1, ...
        phase_init(1), 100, cost_d_up, cost_r_end); save('Data/img1.mat','img1'); 
    img2 = gen_max_freq_arr(n, 0:iterCost:endval, beginH:iterH:1, ...
        phase_init(2), 100, cost_d_up, cost_r_end); save('Data/img2.mat','img2');
    img3 = gen_max_freq_arr(n, 0:iterCost:endval, beginH:iterH:1, ...
        phase_init(3), 100, cost_d_up, cost_r_end); save('Data/img3.mat','img3');
end
cm = brewermap(200,'BuGn');

for i = 1:3
    axes(ha(3+i))
    if i == 1
        img = load('Data/img1.mat','img1'); img = img.img1;
        img(img < 0.01) = 0.01;
    elseif i == 2
        img = load('Data/img2.mat','img2'); img = img.img2;
        img(img < 0.05) = 0.05;
    elseif i == 3
        img = load('Data/img3.mat','img3'); img = img.img3;
        img(img < 0.1) = 0.1;
    end
    imagesc(beginH:iterH:1, 0:iterCost:endval, img);
    set(gca,'clim',[0,1])
    set(gca,'ydir','normal')
    colormap(cm);
    xlabel('Homing efficiency (H)')
    if i == 1
        ylabel('Fitness cost of A')
    elseif i == 3
        pos = get(gca,'position');
        colorbar;
        set(gca,'position',pos)
    end
    set(gca,'xtick',beginH:.1:1)
    set(gca,'ytick',0:.1:endval)
    xlim([beginH, 1])
    ylim([0,endval])
    
    hold on;
    contvals = [0.5, 0.95];
    for cont_i = 1:length(contvals)
        contour(...
            beginH:iterH:1, ...
            0:iterCost:endval,...
            img,...
            [contvals(cont_i), contvals(cont_i)],...
            'w-','linewidth',1 ...
        );
    end
    
    set(gca,plot_opts{:})

end

end

function [arr] = gen_max_freq_arr(n, final_cost_vals, H_vals, init, t_max, d, cost_r)

cl = length(final_cost_vals);
Pl = length(H_vals);

arr = -1*ones(cl,Pl);

ct = 1;
for i = 1:cl
    c = final_cost_vals(i);
    for j = 1:Pl
        H = H_vals(j);
        [~,Y_D,~] = simulate_resistance(n, init, H, c, d, cost_r, t_max, 1);
        A_max = max(Y_D(:,n));
        arr(i,j) = A_max;
    end
    disp(['Done with ' num2str(ct) ' out of ' num2str(cl)])
    ct = ct + 1;
end

disp(['Done with init = ' num2str(init)])

end