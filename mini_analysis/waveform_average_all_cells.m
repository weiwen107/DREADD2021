% This script loops through the waveform_average_kinetics folder and groups
% them by experimental conditions: 
% For example:
    % DR+CNO- 1
    % CNO- 2
    % NT- 3

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  

% Once the waveform averages of all cells have been grouped by their
% experimental conditions, a mega waveform average will be generated by
% taking the average of all cells within each condition.

%% %% Change these accordingly based on how you want to group the data
% in each ALL_EVENTS file, dates can be extracted from the last six digits

%save results
save_results = 1;

%figure on
figure_on_unscaled = 1;
figure_on_scaled = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'WT_CP_mIPSC';

%rise time cutoff
rise = 'mIPSC';

%where to save grouped files
fp_wa_group = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\waveform_average_by_groups\';

%experimental conditions
exp_con = {'DR_CNO','CNO'};

%import cell_id_index table 
cd('C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\mini_data_by_groups')
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%location of waveform average files (per cell)
WA_fp = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\waveform_average_kinetics\';

%% group data by condition

%pre-allocation
all_wavg = cell(1,numel(exp_con));

%loop through folder
cd(strcat(WA_fp, rise))
all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    curr_date = curr_name(1:6);
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(wavg_per_cell{1,1},2)
            if sum(wavg_per_cell{1,1}(:,ci)) == 0
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    all_wavg{1,cond_i}(:,cell_ID) = wavg_per_cell{1,1}(:,ci);
                end
            end
        end
    end
        
end

%% calculate meta average
%pre-allocation
all_meta_wavg = cell(1,numel(exp_con));

%calculate average for each row
for ei = 1:numel(exp_con)
    for ri = 1:size(all_wavg{1,ei},1)
        all_meta_wavg{1,ei}(ri,1) = nanmean(all_wavg{1,ei}(ri,:));
    end
end

%% Scale peaks: scale condition2 to condition1
% usually scale the DR+CNO condition (first cell array) 
% to the CNO condition (second cell array) 

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor = [0 0 0];

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%color code
% color1 = fun([0 0 0],tint_factor(1)); %black-NT
% color2= fun([70 130 180],tint_factor(2)); %steel blue,  tint-24
% color3= fun([70 130 180],tint_factor(3)); %steel blue, 50% tint-48
color2 = fun([250 128 114],tint_factor(2)); %salmon
color1= fun([70 130 180],tint_factor(1)); %steel blue

condition1 = all_meta_wavg{1,1};
condition2 = all_meta_wavg{1,2};
% condition3 = all_meta_wavg{1,1};

[peak_cond1, peak_ind_cond1] = min(condition1);
[peak_cond2, peak_ind_cond2] = min(condition2);

scale_factor = mean(condition1((peak_ind_cond1-1):(peak_ind_cond1+1),1))/...,
    mean(condition2((peak_ind_cond2-1):(peak_ind_cond2+1),1));

condition2_scaled = NaN(size(condition1,1),1);
diff = size(condition1,1)-size(condition2,1);

if diff >=0
    condition2_scaled(diff+1:end) = condition2 .* scale_factor;
else
    condition2_scaled = condition2(abs(diff)+1 : end) .* scale_factor;
end
    

%plotting unscaled

plot_range = (10:235);
y_limit = [-30 5];

if figure_on_unscaled == 1
    figure();
    plot(condition1(plot_range),'Color',color1,'LineWidth',4)
    hold on
    %plot(condition3(plot_range),'Color',color3,'LineWidth',4)   
    plot(condition2(plot_range),'Color',color2,'LineWidth',4)
    %hold on
    %plot(meta_ave.DR_CNO_24(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    %draw scale
    scale_x_start = plot_range(end)-30;
    plot([scale_x_start; scale_x_start+10],[-10; -10], '-k',[scale_x_start;scale_x_start],[-10; -5], '-k', 'LineWidth',2)
    text(scale_x_start, -8, '5 pA', 'HorizontalAlignment', 'right')
    text(scale_x_start+5, -11, '2 ms', 'HorizontalAlignment', 'center')
    
    title('unscaled')
    ylim(y_limit)
    
    box off
    hold off
end

    %plotting scaled
if figure_on_scaled == 1    
    figure();
    plot(condition1(plot_range),'Color',color1,'LineWidth',4)
    hold on
    plot(condition2_scaled(plot_range),'Color',color2,'LineWidth',4)
%     hold on
%     plot(shk3_scaled(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
    title('scaled')
    ylim(y_limit)
    %legend('WT','shk3');
    box off
    hold off
end

%% Fit decay phase of meta averages


%% save to file

if save_results == 1
    cd(strcat(fp_wa_group, rise))
    save(strcat(exp_name,'.mat'),'all_wavg','all_meta_wavg','scale_factor','exp_con')
end
