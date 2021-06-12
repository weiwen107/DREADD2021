%This script records the X coordinates of raw mini data which will be used
%for analysis

%% to be changed for each new run
% where to save excised x coordinates
excised_filename = {'excised_210218.mat'};

% specify experiment
experiment = {'210218'};

% location to saved excised coordinates
fp_excised_pts = 'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\excised_data';

% location of mini data
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\mini\';

%% start excising raw data

cd(strcat(fp_data,experiment{1}))

excised_points = cell(1,numel(experiment)); %numel() returns the number of elements in the array
excised_wavenames = cell(1,numel(experiment)); %cell(n1,n2) returns an n1-by-n2 cell array of empty matrices

for jj = 1:numel(experiment)
    current_experiment = experiment{jj};
    
    all_files = dir; %dir lists files and folders in the current folder
    file_num = numel(all_files)-2; %all data files
    filename = cell(1,file_num);
    cell_id = cell(1,file_num);
    cell_num = NaN(1,file_num);
    filename_bycell = cell(1,file_num);
    
    excised_points{1,jj} = cell(1,file_num);
    excised_wavenames{1,jj} = cell(1,file_num);
    
   %%%%Cell ID index     
   for ff = 1:file_num
            filename{ff} = all_files(ff+2).name;
            num_filename = numel(filename{ff}); %number of characters in the file name

                
            if isnan(str2double(filename{ff}(6))) == 1
                cellID = str2double(filename{ff}(5));
            else
                cellID = str2double(filename{ff}(5:6));
            end

            cell_num(ff) = cellID;
            trace_id = str2double(filename{ff}(num_filename-6 : num_filename-3));


            cell_id{1,cellID}(trace_id,1) = trace_id;
            
            filename_bycell{1,cellID}{trace_id} = filename{ff};
   end
   for ci = 1:max(cell_num)
       if isempty(cell_id{1,ci}) == 1
           continue
       end
       disp(strcat('Current Cell: ', num2str(ci)))
       
       for ti = 1: cell_id{1,ci}(end,1)
           
           disp(strcat('Excising Trace', num2str(ti)))
           trace_id = cell_id{1,ci}(ti,1);
   

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
            if isempty(filename_bycell{1,ci}{ti}) == 1
                disp(strcat('Trace', num2str(ti), ' does not exist!'))
                continue
            else
            
                dtstruct = ws.loadDataFile(filename_bycell{1,ci}{ti});
                prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
                vals = eval(prefix);
                vals = vals(:,1);
            end

            %clear resampled
            clear resampled;
            %resample at 1/2
            resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2); 
            %resample(x,p,q) resamples the input sequence x at p/q times
            %the original sample rate
            %resampled = vals; %remove the downsampling

            figure('position',[119 171 1210 611])
            hold on
            axis([1 numel(resampled) (nanmean(resampled)-120) (nanmean(resampled)+50)])
            %nanmean() returns the sample mean ignoring NaNs
            title(strcat('Cell',num2str(ci), ': trace', sname),'Interpreter', 'none');
            plot(resampled)

            %disp(strcat(current_experiment,'_',num2str(g),' choose start and end points'))
            [X,Y] = ginput; %ginput raises crosshairs in the current axes to identify points in the figure with the mouse
            excised_points{1,jj}{1,ci}{ti} = X;
            excised_wavenames{1,jj}{1,ci}{ti} = filename_bycell{1,ci}{ti};
            close all
            fclose all;
       end %trace
   end %cell
end %experiment


%%
%%%%%%%%%%%%% SAVE  (careful with changes!)

result = input('Save results? This will overwrite previous file unless renamed! (1 = y, 2 = n):','s');

result = str2double(result);

if result == 1
    cd (fp_excised_pts)
    save(excised_filename{1},'excised_points','excised_wavenames')
end
