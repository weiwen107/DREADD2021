% This script extract all amplitudes of mini events from an experiment and 
% save events from each cell separately.

%% Change these before each run
%whether to save results
save_file = 1;

%risetime cutoff
rise = 'rise_1';

%location of mini analysis files
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_mini_results\';

%name of the mini analysis file
mini_name = 'MINIANALYSIS_201029.mat';

%location to save the mini amplitude file
fp_all_events = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\all_mini_events_by_cell\';

%name of the the saved file
save_file_name = strcat('ALL_EVENTS', mini_name(end-10:end));

%% Amplitude extraction

%load file
curr_exp = load(strcat(fp_analyzed_data, rise, '\', mini_name));
curr_exp_amp = curr_exp.AMP_ALL;

cell_num = numel(curr_exp_amp{1,1});

%pre-allocation
all_amp_raw = cell(1,cell_num);
all_amp = cell(1,cell_num);

%extract amplitudes from each trace
tracect = 1;

for ci = 1:cell_num
    for ti = 1:size(curr_exp_amp{1,1}{1,ci},2)
        if sum(~isnan(curr_exp_amp{1,1}{1,ci}(:,ti))) == 0 || ...
                nanmean(curr_exp_amp{1,1}{1,ci}(:,ti)) == 0
            continue
        else       
            lgth = size(curr_exp_amp{1,1}{1,ci}(:,ti),1);
            all_amp_raw{1,ci}(1:lgth,tracect) = curr_exp_amp{1,1}{1,ci}(:,ti);
            tracect = tracect + 1;
        end
    end
end

%convert NaNs into zeros, and then remove zeros before concatnenating all
%columns into a single column

for cii = 1:cell_num
    for row = 1:size(all_amp_raw{1,cii},2)
        for col = 1:size(all_amp_raw{1,cii},1)
            if isnan(all_amp_raw{1,cii}(col,row))
                all_amp_raw{1,cii}(col,row) = 0;
            end
        end
    end
    
    all_amp{1,cii} = nonzeros(all_amp_raw{1,cii});
end

%% save results
if save_file == 1
    cd(strcat(fp_all_events, rise))
    save(save_file_name,'all_amp')
end