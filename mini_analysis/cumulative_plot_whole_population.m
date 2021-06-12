% This script takes grouped mini events, randomly picks 100 events from
% each cell under each experimental condition, and then calculates the
% cumulative distribution

%% load mini populations (grouped by experimental condition)
fp_all_mini_by_group = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\cumulative_data\all_mini_by_group\';
exp_name = 'shk3_KO_CP_DW_24CNO';

rise = 'rise_1';

cd(strcat(fp_all_mini_by_group, rise))
load(strcat(exp_name,'.mat'))

save_results = 1;
fp_cumu = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\cumulative_data\cumulative_stats\';
save_file_name = strcat(exp_name, '_cumu_stats.mat');


%% prepare data
% choose experimental conditions you want to scale and plot for comparison
% usually choose DR+CNO (exp) and CNO (ctrl)(first two cell arrays)
exp_ind = 1; %DR+CNO
ctrl_ind = 2; %CNO

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor = [0 0 0];

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%color code
color1 = fun([70 130 180],tint_factor(1)); %steel blue, 50% tint
color2 = fun([250 128 114],tint_factor(2)); %salmon, 50% tint
%color1 = fun([153 186 221],tint_factor(exp_ind)); %carolina blue
%color2 = fun([255 203 164],tint_factor(ctrl_ind)); %deep peach
%color2 = fun([144 200 144],tint_factor(ctrl_ind)); %light green

exp_amps = selected_mini_events{1,exp_ind}(:,1); 
ctrl_amps = selected_mini_events{1,ctrl_ind}(:,1); 

if size(exp_amps,1) > size(ctrl_amps,1)
    rand_index = randi([1,size(exp_amps,1)],size(ctrl_amps,1),1);
    exp_amps_final = exp_amps(rand_index);
    ctrl_amps_final = ctrl_amps;
    
elseif size(exp_amps,1) < size(ctrl_amps,1)
    rand_index = randi([1,size(ctrl_amps,1)],size(exp_amps,1),1);
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps(rand_index);
    
else
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps;
end    

exp_amps_sort = sort(exp_amps_final);
ctrl_amps_sort = sort(ctrl_amps_final);

%% Scaling

%use linear fitting to get scale factor (scale exp to ctrl)
f = fittype('a*x+b','coefficients',{'a','b'});
slope_start = (max(exp_amps_sort)-min(exp_amps_sort))/(max(ctrl_amps_sort)-min(ctrl_amps_sort));
[lnfit,gof] = fit(ctrl_amps_sort,exp_amps_sort,f,'StartPoint',[slope_start 0]);
exp_scaled = (exp_amps_sort-lnfit.b) ./ lnfit.a;

%plot scaled data
figure()
plot(ctrl_amps_sort,ctrl_amps_sort,'k', 'Linewidth',2)
hold on
plot(ctrl_amps_sort,exp_amps_sort,'o','MarkerSize',8)
plot(lnfit)

ax1 = gca;
ax1.FontSize = 12;
ax1.LineWidth = 2;
ax1.YLabel.String = 'DR+CNO (pA)';
ax1.YLabel.FontSize = 14;
%ax1.XTick = [0 20 40 60 80];
ax1.XLabel.String = 'CNO (pA)';
ax1.XLabel.FontSize = 14;
legend('CNO vs. CNO','DR+CNO vs. CNO','FontSize',12,'Location','southeast')
legend('boxoff')
text(10,max(ctrl_amps_sort),strcat('y=',num2str(lnfit.a),'x',num2str(lnfit.b)),'FontSize',12)
%title('P45-50')
box off

hold off

%% generate cumulative distribution plots

%cumulative distribution for scaled exp data
[cumu_exp_X_scaled, cumu_exp_Y_scaled] = cumhist(exp_scaled,[min(exp_scaled) max(exp_scaled)],0.01);

%two-sample Kolmogorov-Smirnov test
%after scaling
[h1,p1,ks2stat1] = kstest2(exp_scaled,ctrl_amps_sort);
%before scaling
[h2,p2,ks2stat2] = kstest2(exp_amps_sort,ctrl_amps_sort);


figure();
plot(cumulative{1,ctrl_ind}(:,1), cumulative{1,ctrl_ind}(:,2),'Color',color2,'LineWidth',2)
hold on

%plot(cumu_exp_X_scaled, cumu_exp_Y_scaled,'k','LineWidth',2)

plot(cumulative{1,exp_ind}(:,1), cumulative{1,exp_ind}(:,2),'Color',color1,'LineWidth',2)

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'Cumulative %';
ax.YLabel.FontSize = 14;
ax.XLim = [5 45];
ax.YLim = [0 110];
ax.YTick = [0 50 100];
ax.XLabel.String = 'mEPSC amplitude (pA)';
ax.XLabel.FontSize = 14;
legend('CNO','Scaled','DR+CNO','FontSize',12,'Location','southeast')
legend('boxoff')
text(ax.XLim(2)-20,70,strcat('DR+CNO vs.CNO: p=',num2str(p2)),'FontSize',12)
text(ax.XLim(2)-20,60,strcat('Scaled vs.CNO: p=',num2str(p1)),'FontSize',12)
title(strcat(num2str(ax.XLim(2)),'pA cutoff'))
box off
hold off

%% save files
if save_results == 1
    cd(strcat(fp_cumu, rise))
    save(save_file_name,'exp_amps_sort','ctrl_amps_sort','exp_scaled','cumu_exp_X_scaled','cumu_exp_Y_scaled',...
        'h1','p1','ks2stat1','h2','p2','ks2stat2','rand_index')
end
