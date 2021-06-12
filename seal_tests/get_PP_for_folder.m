
%this script is written to be run from a folder of .h5 files
%Passive properties of each trace obtained from each cell will be saved in
%distinct cell arrays and then averages will be calculated

%% Change these for each run
save_results = 1;

experiment = '210603';

fp_pp = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_seal_test';

fp_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\';

% what kind of experiment this is: mode = 1 for mini recordings, mode = 2
% for f-I/sPSC recordings, 3 for sEPSC recordings
experiment_mode = 2;

if experiment_mode == 1
    fp_seal_data = strcat(fp_data, 'mini\',experiment);
elseif experiment_mode == 2
    fp_seal_data = strcat(fp_data, 'seal_test\',experiment);
elseif experiment_mode ==3
    fp_seal_data = strcat(fp_data, 'sp_currents\excitation\',experiment);
end

%parameters of the seal test
step_start = 0.5; %in s

vstep = -0.005; %in V

figure_on = 0; %1 for on, 0 for off

%% looping through the folder

cd(fp_seal_data)
all_files = dir;
file_num = numel(all_files)-2;

sprate = 10000; %rate after resampling

%pre-allocation
filename = cell(1,file_num);
cell_id = cell(1,file_num);
aDAT = cell(1,file_num);
cell_num = NaN(1,file_num);

Rs_est = cell(1,file_num);
Rs = cell(1,file_num);
Rt = cell(1,file_num);
Cm = cell(1,file_num);
Vm = cell(1,file_num);

%Get file names of the specified cell
for f = 1:file_num
    filename{f} = all_files(f+2).name;
    num_filename = numel(filename{f}); %number of characters in the file name
   
    if isnan(str2double(filename{f}(6))) == 1
        cellID = str2double(filename{f}(5));
    else
        cellID = str2double(filename{f}(5:6));
    end
    
    cell_num(f) = cellID;
    trace_id = str2double(filename{f}(num_filename-6 : num_filename-3));
    
    
    cell_id{1,cellID}(trace_id,1) = trace_id;
    
    
    if trace_id < 10
                sname = strcat('000',num2str(trace_id));
            elseif trace_id >= 10 && trace_id < 100
                sname = strcat('00',num2str(trace_id));
            elseif trace_id >= 100 && trace_id < 1000
                sname = strcat('0',num2str(trace_id));
            else
                sname = num2str(trace_id);
    end
    
    
    %datafile readout
    dtstruct = ws.loadDataFile(filename{f});
    prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
    vals = eval(prefix);
    vals = vals(:,1);
    
    %resample @1:2 (sample rate = 5000 Hz)
    %resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);
    
    
    %In aDAT, each cell represents a neuron; within each cell,each column
    %is a single file
    aDAT{1,cellID}(1:size(vals,1),trace_id) = vals;
    
    
end

for gi = 1:max(cell_num)
    cell_id{1,gi} = nonzeros(cell_id{1,gi}); %remove zeros in the cell_id array
    if isempty(cell_id{1,gi}) ==1
        continue
    else
        
        for ri = 1:cell_id{1,gi}(end,1)
            if ismember(ri,cell_id{1,gi}) == 0
                 aDAT{1,gi}(:,ri) = NaN;
             end
         end
    end
end

for ci = 1:max(cell_num)
    if isempty(cell_id{1,ci}) == 1
        continue
    else
        for ti = 1:cell_id{1,ci}(end,1)
            if ismember(ti,cell_id{1,ci}) == 0
                Rs_est{1,ci}(ti,1) = NaN; 
                Rs{1,ci}(ti,1) = NaN;
                Rt{1,ci}(ti,1) = NaN;
                Cm{1,ci}(ti,1) = NaN;
                Vm{1,ci}(ti,1) = NaN;
            else
         
                data = aDAT{1,ci}(:,ti);

                [Rs_est_scaled, Rs_scaled, Rt_scaled, Cm_scaled, Vm_scaled] = ...
                    get_passive_properties(data, step_start, vstep, sprate, figure_on);

                Rs_est{1,ci}(ti,1) = Rs_est_scaled;
                Rs{1,ci}(ti,1) = Rs_scaled;
                Rt{1,ci}(ti,1) = Rt_scaled;
                Cm{1,ci}(ti,1) = Cm_scaled;
                Vm{1,ci}(ti,1) = Vm_scaled;
            end
        end
    end
end

%Calculating averages for each cell
Rs_est_ave = zeros(cellID,1);
Rs_ave = zeros(cellID,1);
Rt_ave = zeros(cellID,1);
Cm_ave = zeros(cellID,1);
Vm_ave = zeros(cellID,1);

for cj = 1:max(cell_num)
    Rs_est_ave(cj,1) = nanmean(Rs_est{1,cj});
    Rs_ave(cj,1) = nanmean(Rs{1,cj});
    Rt_ave(cj,1) = nanmean(Rt{1,cj});
    Cm_ave(cj,1) = nanmean(Cm{1,cj});
    Vm_ave(cj,1) = nanmean(Vm{1,cj});
    
end
    
%% save to file
if save_results == 1 && ~(experiment_mode == 3)
    cd(fp_pp)
    
    display_name = strcat(experiment,'.mat');
    save(display_name,'Rs_est','Rs','Rt','Cm','Vm','Rs_est_ave',...
    'Rs_ave','Rt_ave','Cm_ave','Vm_ave')

end