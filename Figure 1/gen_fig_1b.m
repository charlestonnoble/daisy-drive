function [] = gen_fig_1b()

addpath('../../Plotting Functions/')
addpath('../../')

% Parameters
n       = 3;
d0      = 0.08;
h       = 0.9;
c       = 0.1;
d       = 0.05;
s       = 1;
t_max   = 200;

[T,y_drive,~] = simulate_resistance(n,d0,h,c,d,s,t_max,0);

close all; figure('position',[487   330   298   217],'color','w');
ha = tight_subplot(1,1,[0 0],[0.12,0.05],[0.1,0.1]);
axes(ha); hold on;
cm = brewermap(3,'BrBG');

for i = 1:n
    plot(T,y_drive(:,i),'color',cm(i,:),'linewidth',1)
end
leg = legend('C','B','A'); 
set(leg,'box','off');
ylim([0,1])

end