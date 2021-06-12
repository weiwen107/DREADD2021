%%%% This script takes grouped fI data sets and conducts group analysis

%% Group info
%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\fI_data_by_groups';

%name of the mat file
filename = 'hM4Di_adult_fI_grouped_analysis.mat';

%save results
save_results = 1;

%name of groups (experimental conditions)
cond = {'adult_CNO_24_n_48','adult_DR_CNO','adult_DR_CNO_48h'};

%% import grouped fI data 
%groups are sequentially stored as indicated in cond
group_struct = cell(1,numel(cond));
area_struct = cell(1,numel(cond));

for gi = 1:numel(cond)
    group_struct{1,gi} = eval(cond{1,gi});
    area_struct{1,gi} = eval(strcat(cond{1,gi},'_area'));
end

%% Bursting neurons
% Bursts are defined as two spikes having an interspike interval less than
% 10 ms.

% quantifing IFR at 300pA current injection because: 
% 1. at 200pA, almost all cells have IFRs lower than 100 Hz
% 2. at 400 pA, almost all cells have IFRs higher than 100 Hz
% 300 pA seems like a middle ground
% same reason for choosing 400pA for mean IFR

stim_IFR = 300;
stim_mIFR = 400;

IFR_burst_cell_counter = zeros(1,numel(cond));
IFR_burst_cell_index = cell(1,numel(cond)); %binary indicator: 0 as non-burst, 1 as burst
IFR_prct_burst = zeros(1,numel(cond));

mIFR_burst_cell_counter = zeros(1,numel(cond));
mIFR_burst_cell_index = cell(1,numel(cond)); %binary indicator: 0 as non-burst, 1 as burst
mIFR_prct_burst = zeros(1,numel(cond));

%IFR
for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.IFR,2)
        if isnan(group_struct{1,gii}.IFR(1,ci))
            continue
        else
            curr_IFR = group_struct{1,gii}.IFR(stim_IFR/20,ci);
            if curr_IFR > 100 || curr_IFR == 100
                IFR_burst_cell_counter(1,gii) = IFR_burst_cell_counter(1,gii) + 1;
                IFR_burst_cell_index{1,gii}(1,ci) = 1;
            else
                IFR_burst_cell_index{1,gii}(1,ci) = 0;
            end
        end
    end
    
    IFR_prct_burst(1,gii) = IFR_burst_cell_counter(1,gii) / size(IFR_burst_cell_index{1,gii},2) * 100;
end

%mean IFR
for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.mean_IFR,2)
        if isnan(group_struct{1,gii}.mean_IFR(1,ci))
            continue
        else
            curr_mIFR = group_struct{1,gii}.mean_IFR(stim_mIFR/20,ci);
            if curr_mIFR > 100 || curr_mIFR == 100
                mIFR_burst_cell_counter(1,gii) = mIFR_burst_cell_counter(1,gii) + 1;
                mIFR_burst_cell_index{1,gii}(1,ci) = 1;
            else
                mIFR_burst_cell_index{1,gii}(1,ci) = 0;
            end
        end
    end
    
    mIFR_prct_burst(1,gii) = mIFR_burst_cell_counter(1,gii) / size(mIFR_burst_cell_index{1,gii},2) * 100;
end

%% Correlation between fI curve slope/area under f_I curve and input resistance (using mean IFR)

fI_slope = cell(1,numel(cond));

%r1- correlation coefficients for slope
%p1- corresponding p values for slope
%z1- z transformation of r1

%r2- correlation coefficients for area
%p2- corresponding p values for area
%z2- z transformation of r2

corr_val = cell(1,numel(cond));
z_score = cell(1,numel(cond));

figure_on = 0;

for gii = 1:numel(cond)
    for ci = 1:size(group_struct{1,gii}.mean_IFR,2)
        if isnan(group_struct{1,gii}.mean_IFR(1,ci))
            continue
        else
            curr_mIFR = group_struct{1,gii}.mean_IFR(1:20,ci);
            fit_start = find(curr_mIFR,1,'first');
            curr_X = (fit_start*20:20:400)';
            curr_Y = curr_mIFR(fit_start:20);
            
            curr_fit = fitlm(curr_X,curr_Y);
            fI_slope{1,gii}(ci,1) = curr_fit.Coefficients{2,'Estimate'};
        end
    end
    
    curr_Rin = group_struct{1,gii}.Rin(1:size(fI_slope{1,gii},1));
    
    %slope vs Rin
    [r1,p1] = corrcoef(curr_Rin,fI_slope{1,gii});
    corr_val{1,gii}.r1 = r1;
    corr_val{1,gii}.p1 = p1;
    corr_val{1,gii}.z1 = 0.5*(log(1+r1)-log(1-r1));
    
    %area vs Rin (use mean_IFR area)
    curr_area = area_struct{1,gii}.mean_IFR(1:size(curr_Rin,1));
    [r2,p2] = corrcoef(curr_Rin, curr_area);
    corr_val{1,gii}.r2 = r2;
    corr_val{1,gii}.p2 = p2;
    corr_val{1,gii}.z2 = 0.5*(log(1+r2)-log(1-r2));
    
    
    if figure_on == 1
        
        figure('position',[56 200 800 490])
        subplot(1,2,1)
        scatter(curr_Rin,fI_slope{1,gii})
        hold on
        fitt1 = fitlm(curr_Rin,fI_slope{1,gii});
        plot(fitt1)
        title('Rin vs. fI_slope')
        hold off
        
        subplot(1,2,2)
        scatter(curr_Rin,curr_area)
        hold on
        fitt2 = fitlm(curr_Rin,curr_area);
        plot(fitt2)
        title('Rin vs. fI_area')
        hold off
    end
    
    
end

%% calculate z score for two correlations (2nd-1st)

z_slope(1,1) = (corr_val{1,2}.z1(2,1)-corr_val{1,1}.z1(2,1))/...,
          sqrt(1/(size(group_struct{1,2}.MFR,2)-3) +...,
          1/(size(group_struct{1,1}.MFR,2)-3));
%calculate p value
z_slope(1,2) = (1-normcdf(z_slope(1,1)))*2;

if abs(z_slope(1,1))>1.96 %critical z value for 95% confidence interval
    z_slope(1,3) = 1; %reject null hypothesis, two correlation differnt
else
    z_slope(1,3) = 0; %no siginificant change in correlation robustness
end

z_area(1,1) = (corr_val{1,2}.z2(2,1)-corr_val{1,1}.z2(2,1))/...,
          sqrt(1/(size(group_struct{1,2}.MFR,2)-3) +...,
          1/(size(group_struct{1,1}.MFR,2)-3));
      
z_area(1,2) = (1-normcdf(z_area(1,1)))*2;

if abs(z_area(1,1))>1.96 %critical z value for 95% confidence interval
    z_area(1,3) = 1; %reject null hypothesis, two correlation differnt
else
    z_area(1,3) = 0; %no siginificant change in correlation robustness
end
      
%% save results

if save_results == 1

    cd (fp_analyzed_data)
    save(filename,'stim_IFR','stim_mIFR','IFR_burst_cell_counter','IFR_burst_cell_index','IFR_prct_burst',...
        'mIFR_burst_cell_counter','mIFR_burst_cell_index','mIFR_prct_burst',...
        'fI_slope','corr_val','z_slope','z_area','cond')
end

