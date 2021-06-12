%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% mEPSC ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What you will need to run this script:

% 1) Organize files.  Cell data should be organized in a folder structure
% on your computer by experiment. Files belonging to a specific cell
% should have the same prefix "cellx' where x indicates the cell number.
% Files associated with this script as well as this script
% should be kept in a mini analysis folder within your MATLAB folder.

% 2) The script automatically generate a cell_id cell array which serves as
% a file management system. Each file should be named as such:
% cell1_0001.h5. The cell_id array stores the experiment number, the cell
% number, and the trace number hierarchically, for example, cell_id{1,1}
% indicates experiment 1, under which you will find several cell number arrays
% representing each cell (cell_id{1,1}{1,1} indicates the cell1 in
% experiment 1). Within each cell number array, the id # of all traces that belongs to
% this cell is stored in a column.

% 3) Excise X data.  Run the included script excise_x_points* to manually
% remove seal tests, noisy regions, and noisy traces.  This only needs to
% be done once.

%% these need to be changed for each new data run

%Name of file you would like to have data saved to
save_file_name = 'MINIANALYSIS_201026.mat';

%this needs to be copied from excise_x_points
excised_filename = {'excised_201026.mat'};

%loop through all the folders
experiment = {'201026'};

%location of the excised file
excised_data_fp = strcat('C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\excised_data\',excised_filename{1});

%Filepath where you keep data folders
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\chronic_hm4di\mini\';

%location to save the analyzed mini data
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_mini_results\';

%% parameters for analysis

figures_on = 1;
%%%%%% 0 = NO, 1 = YES

troubleshoot_events_detect = 0;
%%%%%% opens a figure window that shows individual events that have been
%%%%%% selected/analyzed and displays relevant event properties

troubleshoot_events_wavg = 0;
%%%%%% opens a figure window that shows individual events selected for wavg
%%%%%% analysis and decay fitting

%%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

samplert = 5000;
%%%%%% Sampling rate at which to analyze data

tRISE = 0.0005*samplert;

tDECAY = 0.0025*samplert;

num_zeros_temp_BL = 10;

num_pts_temp = 18;
%%%%%% Number of points in template

detect_crit_thresh = 0.10;
%%%%%% Threshold value for 'acceptable' template match w/ data

%%%%%%%%%%%%%%%%%%%%%%%%%% EXCLUSION CRITERIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

risetime_cutoff = 1;
%%%%%% (ms) Detected events must have SHORTER OR EQUAL rise time than this value (in order to be accepted)
%%%%%% NOTE: It's easy to remove events post-analysis with high rise times,
%%%%%% so probably best to keep this high (3+)

amp_cutoff_low = 5;
amp_cutoff_high = 120;
%%%%%% (pA) Value range of acceptable events amplitudes

refractory_per = 15;%15;
%%%%%% minumum number of points between events (REFRACTORY PERIOD)

numptsBL = 10;
%%%%%% number of points before event start to use for determining whether
%%%%%% event has a stable baseline

BL_slope_thresh_hi = 0.8;
BL_slope_thresh_lo = -0.8;
%%%%%% Slope of baseline period before individual events must be between bounds

decay_frac_thresh = 0.9;
%%%%%% end of decay is defined as point at which smoothed data returns to
%%%%%% this fraction of baseline amplitude

event_charge_thresh = 1e-2;
%%%%%% individual events must have charge (area under curve) greater than
%%%%%% this value to be accepted

fit_noise_thresh = 3.0;%1.9;
%%%%%% If max value of eventBL >= fit_noise_thresh * std(data), exclude event
%%%%%% prevents fitting too many noise transients

BLBAD_thresh = 3.0;
%%%%%% If event start value is more than this much over baseline, exclude
%%%%%% event

num_pts_decay_fit = 30;%35;
%%%%%% number of points for exponential fit to decay

min_decay_rsquare = 0.25;
%%%%%% Minimum R-squared value acceptable for exponential fit to event
%%%%%% decay 

max_decay_rmse = 3;%2.5;
%%%%%% Maximum RMSE value acceptable for exponential fit to event
%%%%%% decay (exclude if higher)

S_E_L = 150;
%%%%%% Standard event length, MUST BE DIVISIBLE BY 3

%% Begin analysis

% -------------------------------------
% TEMPLATE
% NORM is the normalization factor
clear template

NORM            = 0.02;
spoints         = 1: (samplert * 0.020);
template        = zeros(1,samplert * 0.020);

tcount = 0;

for t = 0:length(spoints)
    tcount = tcount+1;
    template(tcount)     = NORM *(1 - exp(-t/tRISE)) * exp(-t/tDECAY); %create a normalized template of an event
end

template = [zeros(1,num_zeros_temp_BL) template]; %%% this forces a stable baseline prior to each mini
template = template(1,1:num_pts_temp);

[~,pk_ind] = max(template);   %%%peak index of template

temp_time = (0:0.2:50);%time points also normalized
temp_time = temp_time(1,1:length(template));
%             figure
%             plot(temp_time,template)
%             axis([0 5 0 max(template)])
n = length(template);
% END TEMPLATE
% -------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIPASS FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Excised Data is run through a high pass filter to remove slow
%%%%%%% deviations from baseline

hfilt = designfilt('highpassiir', ...       % Response type
    'StopbandFrequency',0.35, ...     % Frequency constraints
    'PassbandFrequency',0.75, ...
    'PassbandRipple',0.5, ...
    'StopbandAttenuation',3, ...    % Magnitude constraints
    'SampleRate',5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excised_data = load(excised_data_fp);
excised_points = excised_data.excised_points;
numexp = numel(experiment); 

AMP_ALL = cell(1,numexp);
WAVG = cell(1,numexp); 
TOTAL_TIME = cell(1,numexp);
RISE = cell(1,numexp);
DECAY = cell(1,numexp);
wavgs_raw = cell(1,numexp);
TIME_INDICES = cell(1,numexp);

aDAT_all = cell(1,numexp);
EVENT_START_INDICES = cell(1,numexp);
PT3_AVERAGES = cell(1,numexp);
SMOOTH_PEAK_INDICES = cell(1,numexp);
DECAY_INDICES = cell(1,numexp);
PT3_DECAY_VALUES = cell(1,numexp);

%summary of averages for each mini file in each experiment
meanAmp = cell(1,numexp);
mini_count = cell(1,numexp);
total_time = cell(1,numexp);
FRQ = cell(1,numexp);


for jj = 1:1%numel(experiment) %loop through all folders
    current_experiment = experiment{jj};
    current_folder = strcat(fp_data,current_experiment);
    cd(current_folder)
    
    all_files = dir;
    file_num = numel(all_files)-2;
    filename = cell(1,file_num);
    cell_id = cell(1,file_num);
    cell_num = NaN(1,file_num);
    filename_bycell = cell(1,file_num);
    
    
   %%%% Cell ID index
    for hh = 1:file_num
        filename{hh} = all_files(hh+2).name;
        num_filename = numel(filename{hh}); %number of characters in the file name


        if isnan(str2double(filename{hh}(6))) == 1
            cellID = str2double(filename{hh}(5));
        else
            cellID = str2double(filename{hh}(5:6));
        end
    
        cell_num(hh) = cellID;
        trace_id = str2double(filename{hh}(num_filename-6 : num_filename-3));


        cell_id{1,cellID}(trace_id,1) = trace_id;
        
        filename_bycell{1,cellID}{trace_id} = filename{hh};
    end
    %% loop through all cells
    for ci = 7:7%1 : max(cell_num)
        if isempty(cell_id{1,ci}) == 1
            continue
        else
        
            disp(strcat('Analyzing Cell', num2str(ci)))

            for kk = 1:cell_id{1,ci}(end,1)
                disp(strcat('Analyzing Trace', num2str(kk)))      
                trace_id = cell_id{1,ci}(kk,1);


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
               if isempty(filename_bycell{1,ci}{kk}) == 1
                   disp(strcat('Trace',num2str(kk), ' does not exist!'))

                   FRQ{1,jj}{1,ci}(kk,1) = NaN;
                   AMP_ALL{1,jj}{1,ci}(:,kk) = NaN;

                   meanAmp{1,jj}{1,ci}(kk,1) = NaN;
                   mini_count{1,jj}{1,ci}(kk,1) = NaN;

                   TOTAL_TIME{1,jj}{1,ci}(kk,1) = NaN;
                   RISE{1,jj}{1,ci}(:,kk) = NaN;
                   DECAY{1,jj}{1,ci}(:,kk) = NaN;
                   wavgs_raw{1,jj}{1,ci}{kk} = NaN;
                   TIME_INDICES{1,jj}{1,ci}(:,kk) = NaN;
                   aDAT_all{1,jj}{1,ci}(:,kk) = NaN;
                   continue

               else

                dtstruct = ws.loadDataFile(filename_bycell{1,ci}{kk});
                prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
                vals = eval(prefix);
                vals = vals(:,1);



                %check these
                clear resampled
                clear yy
                clear aDAT
                clear timeindx
                clear n_ind


                samprate = 10000;
                if samprate == 5000  %%%%% if sampling rate = 5000, resample @1:1
                    resampled(1:numel(resample(vals,1,1)),1) = resample(vals,1,1);
                elseif samprate == 10000  %%%%% if sampling rate = 10000, resample @1:2
                    resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);
                end


                yy = [];
                start_end = excised_points{1,jj}{1,ci}{kk}; 
                start_end = start_end(isnan(start_end)==0);
                start_end = start_end(start_end>0);
                if mod(numel(start_end),2) == 0
                    for m = 1:2:numel(start_end-1)
                        if start_end(m+1) > numel(resampled)
                            start_end(m+1) = numel(resampled);
                        end
                        warning('off','all')  %%%%% for some reason, some selected excision points are not integers...
                        %%%%%%%% FILTER DATA
                        yy(numel(yy)+1:numel(yy)+numel(resampled(start_end(m):start_end(m+1)))) = ...
                            filter(hfilt,detrend(resampled(start_end(m):start_end(m+1))));
                        %plot(filter(hfilt,detrend(resampled(start_end(ll):start_end(ll+1)))))
                        warning('on','all')
                    end
                else
                    disp(strcat('Warning...',current_experiment,'_','trace',sname,...
                        ' does not contain an appropriate number of excision points'))
                    yy = filter(hfilt,detrend(resampled));
                end

                yy = yy(isfinite(yy));
                aDAT(1:numel(yy),1) = yy; %excised and filtered data stored in a column (aDAT)
            
                %%%%%%%%%%%%%%%%%% TEMPLATE MATCHING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Cycle through the data and compare the template to each sample point +
                % length template
                clear SCALE OFFSET SSE STANDARDERROR DETECT_CRIT DETECT_CRIT2
                SCALE = NaN(1,numel(aDAT));
                OFFSET = NaN(1,numel(aDAT));
                SSE = NaN(1,numel(aDAT));
                STANDARDERROR = NaN(1,numel(aDAT));
                DETECT_CRIT = NaN(1,numel(aDAT));

                for ii = 1:length(aDAT) - n % difference between data and template

                    dat                 = aDAT(ii:ii+n-1); % resampled and filtered
                    sTop                = sum(template*dat) - sum(template)*sum(dat)/n;
                    sBottom             = sum(template.^2) - sum(template) * sum(template)/n;
                    SCALE(ii)           = sTop/sBottom;

                    OFFSET(ii)          = (sum(dat) - SCALE(ii)*sum(template))/n;

                    if SCALE(ii) > 0
                        fitted_template = zeros(1,numel(template));
                    else
                        fitted_template     = template*SCALE(ii) + OFFSET(ii);
                    end
                     
                    %amp cutoff for templates are 5 and 200
                    if abs(fitted_template(1,1)-fitted_template(1,pk_ind)) > 5 && abs(fitted_template(1,1)-fitted_template(1,pk_ind)) < 200
                        SSE(ii)             = sum(dat)^2 + SCALE(ii)^2 * sum(template.^2) + n*OFFSET(ii)^2 - ...
                            2*(SCALE(ii) * sum(template*dat)) + OFFSET(ii) * sum(dat) - SCALE(ii) * OFFSET(ii) * sum(template);
                        STANDARDERROR(ii)   = (SSE(ii)/(n-1))^(1/2);
                        DETECT_CRIT(ii)     = SCALE(ii)/STANDARDERROR(ii);
                    else
                        SSE(ii)             = sum(dat)^2 + SCALE(ii)^2 * sum(template.^2) + n*OFFSET(ii)^2 - ...
                            2*(SCALE(ii) * sum(template*dat)) + OFFSET(ii) * sum(dat) - SCALE(ii) * OFFSET(ii) * sum(template);
                        STANDARDERROR(ii)   = (SSE(ii)/(n-1))^(1/2);
                        DETECT_CRIT(ii)     = 0;
                    end
                end

                DETECT_CRIT2 = DETECT_CRIT;
                DETECT_CRIT2(DETECT_CRIT2>0) = 0;
                DETECT_CRIT2 = DETECT_CRIT2/min(DETECT_CRIT2);

                thresh          = detect_crit_thresh;     %THRESHOLD FOR DETECTION
                past            = DETECT_CRIT2>thresh;
                hits            = 0.*DETECT_CRIT2;
                hits(past)      = 1;
                hits            = [ 0 diff(hits) ];
                hits_temp = hits;
                hits(hits<0)    = 0;    % finds all points where there is a crossing

                timeindx        = find(hits); %detected events index

                % remove NaNs in aDAT               
                n_ind = find(~isnan(aDAT));
                timeindx = timeindx(ismember(timeindx,n_ind));


                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%% Eliminate false positive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                disp('Detecting events...')

                % remove events detected in the last 350 points to prevent
                % error in subsequent detections
                timeindx(timeindx >= max(timeindx)-200) = NaN; 
                timeindx = timeindx(~isnan(timeindx));
                
                timeindx_ampcutoff = NaN(1,numel(timeindx));
                timeindx_slopeBL = NaN(1,numel(timeindx));
                timeindx_rise = NaN(1,numel(timeindx));
                timeindx_decay = NaN(1,numel(timeindx));
                timeindx_decayfit = NaN(1,numel(timeindx));
                timeindx_charge = NaN(1,numel(timeindx));
                timeindx_BLBAD = NaN(1,numel(timeindx));
                timeindx_noiseSTD = NaN(1,numel(timeindx));
                timeindx_multipk = NaN(1,numel(timeindx));
                timeindx_PeakAtEnd = NaN(1,numel(timeindx));
                timeindx_StartAfterPeak = NaN(1,numel(timeindx));
                timeindx_otherC = NaN(1,numel(timeindx));
                timeindx_LatePeak = NaN(1,numel(timeindx));
                timeindx_refractory_per = NaN(1,numel(timeindx));
                timeindx_decay_start = NaN(1,numel(timeindx));
                timeindx_BLnopts = NaN(1,numel(timeindx));
                timeindx_no_decay = NaN(1,numel(timeindx));

                mean_slope_BL = NaN(1,numel(timeindx));
                event = NaN(S_E_L+1,numel(timeindx));
                norm_event = NaN(S_E_L+1,numel(timeindx));
                pt3_avg = NaN((S_E_L/3)+1,numel(timeindx));
                slope = NaN(S_E_L/3,numel(timeindx));
                slope3_avg = NaN((S_E_L/3)-2,numel(timeindx));
                sm_pk_val = NaN(1,numel(timeindx));
                sm_pk_ind = NaN(1,numel(timeindx));
                slope_pk_val = NaN(1,numel(timeindx));
                slope_pk_ind = NaN(1,numel(timeindx));
                event_start_ind = cell(1,numel(timeindx));
                pt3_decay = NaN(S_E_L*(2/3),numel(timeindx));
                decay_ind = NaN(1,numel(timeindx));%decay start point
                TOTAL_RISE_TIME = NaN(1,numel(timeindx));%from start to peak
                RISE_TIME = NaN(1,numel(timeindx)); %10-90
                DECAY_TIME = NaN(1,numel(timeindx));
                RAW_AMP = NaN(1,numel(timeindx));

                aaDAT = aDAT;

                %%%Fitting function of the decay phase (single exp)
                s = fitoptions('Method','NonlinearLeastSquares','Startpoint',[-11 -0.4],...
                            'Lower',[-400,-3],'Upper',[-5,-0.04]);
                            f = fittype('a*exp(b*x)','options',s);

                for i = 1:numel(timeindx)
                    event(:,i) = aDAT(timeindx(i):timeindx(i)+S_E_L);
                    for k = 4:((S_E_L/3)+1)
                        pt3_avg(k-2,i) = (event(k-2,i)+event(k-3,i)+event(k-1,i))/3;%average amplitudes of every 3 points
                    end
                    pt3_avg(1,i) = (event(1,i)+event(2,i))/2;%first two points
                    pt3_avg(S_E_L/3,i) = (event(S_E_L+1,i)+event(S_E_L,i)+event(S_E_L-1,i))/3;
                    pt3_avg(S_E_L/3+1,i) = (event(S_E_L+1,i)+event(S_E_L,i))/2;
                    for k = 1:S_E_L/3
                        slope(k,i) = pt3_avg(k,i) - pt3_avg(k+1,i); %slope of every 6 points (pre-after)
                    end
                    for k = 4:((S_E_L/3)+1)
                        slope3_avg(k-2,i) = (slope(k-2,i)+slope(k-3,i)+slope(k-1,i))/3;  %averaging slopes
                    end
                    slope3_avg(1,i) = (slope(1,i)+slope(2,i))/2;
                    slope3_avg((S_E_L/3)-1,i) = (slope((S_E_L/3)-1,i)+slope(S_E_L/3,i))/2;

                    %%%find the peak amplitude (value and location)
                    [sm_pk_val(i),sm_pk_ind(i)] = min(pt3_avg((1:end-10),i));
                    %raw_sm_pk_val(i) = sm_pk_val(i);
                    %raw_sm_pk_ind(i) = sm_pk_ind(i);


                    %%%%%%%%%%%%%%%%%%% ELIMINATING FALSE-POSITIVES %%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if sm_pk_val(i) == pt3_avg(end,i)
                        %In case peak is detected at the very end
                        sm_pk_val(i) = NaN;
                        sm_pk_ind(i) = NaN;
                        timeindx_PeakAtEnd(i) = timeindx(i);
                        continue
                    else
                        %%%%%  Event start defined as the point where slope
                        %%%%%  is 30% of the maxium or the first point
                        %%%%%  (when can't find one)
                        [slope_pk_val(i),slope_pk_ind(i)] = max(slope3_avg(1:(end-10),i));
                        event_start_ind{i} = find(slope3_avg(1:slope_pk_ind(i),i) <= 0.3* slope_pk_val(i),1,'last');   
                        if isempty(event_start_ind{i}) == 1
                            event_start_ind{i} = 1;
                        end
                        
                        %reject if event starts after peak
                        if isnan(sm_pk_ind(i)) == 0 && event_start_ind{i} > sm_pk_ind(i)
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_StartAfterPeak(i) = timeindx(i);
                            continue
                        end

                        %reject if peak shows up too late
                        if isnan(sm_pk_ind(i)) == 0 && isempty(event_start_ind{i})==0
                            if event_start_ind{i} >= numptsBL+20
                                sm_pk_val(i) = NaN;
                                sm_pk_ind(i) = NaN;
                                timeindx_LatePeak(i) = timeindx(i);
                            end
                        end
                        
                        %reject multipeak events
                        if isnan(sm_pk_ind(i)) == 0
                            if max(slope3_avg(:,i)) <= 1.5*max(slope3_avg(1:event_start_ind{i}-1,i))
                                sm_pk_val(i) = NaN;
                                sm_pk_ind(i) = NaN;
                                timeindx_multipk(i) = timeindx(i);
                                continue
                            end
                        end  
                    end
                    
                    %the location of the event start point in raw data
                    event_start_aDAT_ind = timeindx(i)+event_start_ind{i}-1;
                    if isnan(sm_pk_ind(i)) == 0
                        if (event_start_aDAT_ind-numptsBL+1) > 0
                            BL_vals = aDAT(event_start_aDAT_ind-numptsBL+1:event_start_aDAT_ind);
                        else
                            %reject if not enough points in baseline
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_BLnopts(i) = timeindx(i);
                        end
                    end

                    %reject if noise exceeds upper limit
                    if isnan(sm_pk_ind(i)) == 0 && ...
                            max(BL_vals-(nanmean(BL_vals(1:numptsBL-3)))) >= fit_noise_thresh*std(aDAT)+mean(aDAT)
                        sm_pk_val(i) = NaN;
                        sm_pk_ind(i) = NaN;
                        timeindx_noiseSTD(i) = timeindx(i);
                        continue              
                    end

                    if isnan(sm_pk_ind(i)) == 0
                        
                        %rejct if event start deviates too much from the
                        %baseline
                        if pt3_avg(event_start_ind{i},i) - nanmean(BL_vals) >= BLBAD_thresh
                            %disp(pt3_avg(event_start_ind{i},i));
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_BLBAD(i) = timeindx(i);
                            continue
                        end
                        
                        if isnan(sm_pk_ind(i)) == 0
                            BL_fit = fitlm((1:numptsBL)',BL_vals);
                            intercept = BL_fit.Coefficients{1,1};
                            slope_fit = BL_fit.Coefficients{2,1};
                            x_coords_fit = ([1 numptsBL]);
                            y_coords_fit = ([slope_fit+intercept slope_fit*numptsBL+intercept]);
                            
                            if troubleshoot_events_detect == 1
                                plot(x_coords_fit,y_coords_fit,'r-');
                            end
                            
                            rsquared = BL_fit.Rsquared.Ordinary;
                            mean_slope_BL(i) = slope_fit; %baseline slope

                            %reject if baseline slope out of bounds
                            if mean_slope_BL(i) < BL_slope_thresh_hi && mean_slope_BL(i) > BL_slope_thresh_lo
                                if troubleshoot_events_detect == 1
                                    shg
                                    plot(x_coords_fit,y_coords_fit-nanmean(BL_vals),'r-');
                                    disp(strcat('event',num2str(i)))
                                    disp(strcat('mean slope BL =',num2str(mean_slope_BL(i))))
                                    disp(strcat('rsquared =',num2str(rsquared)))
                                end
                            else
                                sm_pk_val(i) = NaN;
                                sm_pk_ind(i) = NaN;
                                timeindx_slopeBL(i) = timeindx(i); 
                                continue
                            end
                        end
                    end


                    %%%%% RISE MEASUREMENT
                    if isnan(sm_pk_ind(i)) == 0
                        amplitude = -1*(pt3_avg(sm_pk_ind(i),i) - nanmean(BL_vals));
                        peak_i = sm_pk_ind(i)-event_start_ind{i}+numptsBL;
                        event_wBL = aDAT(event_start_aDAT_ind-numptsBL+1:event_start_aDAT_ind+S_E_L+1-numptsBL);
                        %%%%% event with BL length specified 
                        norm_event(:,i) = event_wBL - nanmean(event_wBL(1:numptsBL));
                        %%%%% normalized to BL avg
                        rise_start_10pc_i = find(norm_event(11:peak_i,i) <= -0.1*amplitude,1,'first');
                        %%%%% 10% rise is first point greater than 10% peak
                        rise_end_90pc_i = find(norm_event(11:peak_i,i) >= -0.9*amplitude,1,'last');
                        %%%%% 90% rise is last point less than 90% peak
                        
                        %reject if cannot find rise range
                        if isempty(rise_start_10pc_i)||isempty(rise_end_90pc_i) == 1
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_rise(i) = timeindx(i);
                            RISE_TIME(i) = NaN;
                        else
                        
                            %%%%% 10-90% rise time region
                            rise_vals_10to90 = norm_event(rise_start_10pc_i:rise_end_90pc_i,i);
                            
                            %0-100% rise time
                            TOTAL_RISE_TIME(i) = (sm_pk_ind(i)-event_start_ind{i})/(samplert);

                            %%%%% 10-90% rise time
                            RISE_TIME(i) = (rise_end_90pc_i-rise_start_10pc_i)/(samplert);
                           
                            if troubleshoot_events_detect == 1
                                disp(strcat('RISE = ',num2str(RISE_TIME(i))))
                            end
                            
                            %reject if rise time exceeds upper limit
                            if RISE_TIME(i) > (risetime_cutoff/1000)
                                sm_pk_val(i) = NaN;
                                sm_pk_ind(i) = NaN;
                                timeindx_rise(i) = timeindx(i);
                                RISE_TIME(i) = NaN;
                                continue
                            end
                        end
                    end


                    %%%%% DECAY MEASUREMENT
                    if isnan(sm_pk_ind(i)) == 0
                        decay_amp = NaN(1,(S_E_L*(2/3)));
                        
                        for k = 2:(S_E_L*(2/3)+1)
                            pt3_decay(k-1,i) = (event(k-1+sm_pk_ind(i),i)+...,
                                event(k+sm_pk_ind(i),i)+event(k+1+sm_pk_ind(i),i))/3;
                            decay_amp(k-1) = nanmean(BL_vals) - pt3_decay(k-1,i);
                        end
                        
                        decay_start_i = find(norm_event(peak_i:peak_i+20,i) >= -0.9*amplitude,1,'first')-1+peak_i;
                        
                        %reject if cannot find decay start point
                        if isempty(decay_start_i) == 1
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_decay_start(i) = timeindx(i);
                            continue
                        end
                        
                        decay_vals_forfit = norm_event(decay_start_i:decay_start_i+num_pts_decay_fit,i);
                        
                        %decay end point
                        decay_i = find(decay_amp <= (1-decay_frac_thresh)*(-1*(pt3_avg(sm_pk_ind(i),i) - nanmean(BL_vals))),1,'first');
                        decay_i_wBL = sm_pk_ind(i)-event_start_ind{i}+numptsBL+2+decay_i;

                        %%%%%% reject if decay end point cannot be found
                        if isnan(sm_pk_ind(i)) == 0 && isempty(decay_i) ==1
                            decay_i = NaN;
                            DECAY_TIME(i) = NaN;
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_no_decay(i) = timeindx(i);
                            
                        else
                            DECAY_TIME(i) = (decay_i+2)/samplert;
                        end
                        
                        %reject if rise is longer than decay
                        if isnan(DECAY_TIME(i)) == 0 && TOTAL_RISE_TIME(i) >= DECAY_TIME(i)
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_decay(i) = timeindx(i);
                            DECAY_TIME(i) = NaN;
                            continue
                        end
                    else
                        decay_i = NaN;
                    end
                    
                    %store the decay end point for each event
                    decay_ind(i) = decay_i;

                    %%%%% FIT DECAY TO EXPONENTIAL
                    if isnan(sm_pk_ind(i)) == 0 
                        decay_vals_forfit = smooth(decay_vals_forfit);
                        [exp_fit,gof] = fit((0:0.2:0.2*num_pts_decay_fit)',decay_vals_forfit,f);

        %                 if timeindx(i) == 17921
        %                 figure
        %                 plot(0:0.2:0.2*num_pts_decay_fit,decay_vals_forfit)
        %                 hold on
        %                 plot(exp_fit)
        %                 disp(gof)
        %                 disp(decay_start_i)
        %                 sm_pk_val(i)
        %                 hold off
        %                 figure
        %                 plot(1:151,norm_event(:,i))
        %                 end

                        %reject if decay is not a good exponential fit
                        if gof.rsquare < min_decay_rsquare || gof.rmse > max_decay_rmse ...
                                || exp_fit.b > -0.03000001
%                             gof.rsquare
%                             gof.rmse

                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            decay_vals_forfit = NaN; 
                            timeindx_decayfit(i) = timeindx(i);
                            continue
                        end
        %                         figure
        %                         hold on
        %                         plot((decay_start_i:decay_start_i+num_pts_decay_fit)',decay_vals_forfit)
        %                         plot(exp_fit)
                    end




                    %%%%% CHARGE MEASUREMENT
                    if isnan(sm_pk_ind(i)) == 0
                        if isnan(decay_i) == 0
                            event_charge = sum(norm_event(numptsBL+1:decay_i_wBL,i))/samplert;
                        else
                            event_charge = sum(norm_event(numptsBL+1:numptsBL+30,i))/samplert;
                        end
                        
                        %reject if charge less than threshold value
                        if -event_charge < event_charge_thresh
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_charge(i) = timeindx(i);
                            continue
                        end
                    end
                    
                   %%%%%%% refractory period between events 
                   
                   % reject if refractory period too low
                   if isnan(sm_pk_ind(i)) == 0  
                       if i>1
                           events_before_i = sm_pk_ind(1:i-1);
                           i_notnan = find(isnan(events_before_i)==0,1,'last'); %distance from the last eligible event
                        if (timeindx(i)+event_start_ind{i})...,
                                - (timeindx(i_notnan)+sm_pk_ind(i_notnan)+decay_ind(i_notnan)) <= refractory_per
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            timeindx_refractory_per(i) = timeindx(i);
                            continue
                        end              
                       end
                   end
                   
                   %trobleshoot individual events
                    if troubleshoot_events_detect == 1
                        hold on
                        plot(1:151,norm_event(:,i))
                        plot(11,0,'go')
                        plot(rise_start_10pc_i:rise_end_90pc_i,rise_vals_10to90,'g-')
                        plot(decay_start_i:decay_start_i+num_pts_decay_fit,decay_vals_forfit,'r-')
                        plot(exp_fit,'y-')

                        if isnan(decay_i) == 0
                            plot(decay_i_wBL,norm_event(decay_i_wBL,i),'ro')
                        else
                            plot(numptsBL+30,norm_event(numptsBL+30,i),'ro')
                        end
                        if isnan(decay_i) == 0
                            for g = numptsBL+1:decay_i_wBL
                                plot([g g],[0 norm_event(g,i)],'y-')
                            end
                        else
                            for g = numptsBL+1:numptsBL+30
                                plot([g g],[0 norm_event(g,i)],'y-')
                            end
                        end
                        disp(strcat('charge =',num2str(-event_charge)))
                    end


                    %%%%% EVENT AMPLITUDE MEASUREMENT
                    if isnan(sm_pk_ind(i)) == 0
                        RAW_AMP(i) = pt3_avg(sm_pk_ind(i),i) - nanmean(BL_vals);
                        
                        %reject if amp out of bounds
                        if -1*RAW_AMP(i) < amp_cutoff_low || -1*RAW_AMP(i) > amp_cutoff_high
                            sm_pk_val(i) = NaN;
                            sm_pk_ind(i) = NaN;
                            RAW_AMP(i) = NaN;
                            timeindx_ampcutoff(i) = timeindx(i);
                        end
                        
                        if troubleshoot_events_detect == 1
                            disp(strcat('amp =',num2str(RAW_AMP(i))))
                            plot(sm_pk_ind(i)-event_start_ind{i}+numptsBL+1,RAW_AMP(i),'bo')
                            hold off
                            waitforbuttonpress
                        end
                    end

                end

                RAW_AMP = -1*RAW_AMP;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAVEFORM AVG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Calculating waveform average...')

                wavgs = NaN(S_E_L,numel(timeindx));
                
                for ii = 1:numel(timeindx)
                    if isfinite(RAW_AMP(ii))==1
                        wave = norm_event(:,ii);
                        wave = wave(1:(S_E_L-20)); %length of each wave used for smoothing
                        sm_wave = smooth(wave);
                        
                        time = (0:0.2:0.2*numel(wave(sm_pk_ind(ii)-event_start_ind{ii}+numptsBL+1:numel(wave))))';
                        y = -sm_wave(sm_pk_ind(ii)-event_start_ind{ii}+numptsBL-3:numel(sm_wave)-3);
                        y = y-mean(y(20:numel(y))); 
                        
                        if y(2,1)>y(1,1)
                            y = y(2:numel(y),1);
                            time = time(2:numel(time),1);
                        end

                        if ii<numel(timeindx)
                            wave_event_start_ind = numptsBL+1;
                            wave_event_pk_ind = sm_pk_ind(ii)-(event_start_ind{ii}-numptsBL);
                            %sm_pk_ind(ii)-event_start_ind{i}+numptsBL+1;
                            rise_phase_ind = find(sm_wave==max(sm_wave(wave_event_start_ind:wave_event_pk_ind)))...
                                :find(sm_wave==min(sm_wave(wave_event_start_ind:wave_event_pk_ind)));
                            if isempty(rise_phase_ind) == 0
                                %     find((sign(sm_wave(rise_phase_ind)-sm_pk_val(ii)/2))==-1,1,'first')
                                tmp_midrise = abs(sm_wave(rise_phase_ind)-sm_pk_val(ii)/2);
                                %     tmp_midrise = abs(sm_wave((sm_pk_ind(ii)-4):sm_pk_ind(ii))-sm_pk_val(ii)/2);
                                [~, mid_rise_ind] = min(tmp_midrise);
                                mid_rise_ind = mid_rise_ind+rise_phase_ind(1,1)-1;
                                wavgs((num_zeros_temp_BL*4)-mid_rise_ind:(numel(wave)-mid_rise_ind+((num_zeros_temp_BL*4)-1)),ii) = wave;
                               
                                if mid_rise_ind >= (num_zeros_temp_BL*4)-1
                                    continue
                                end
                                
                                if troubleshoot_events_wavg == 1
                                    plot(time,y)
                                    hold on
                                    %plot(time,exp_fit(y),'g-')
                                    hold off
                                    waitforbuttonpress;
                                end
                            else
                                disp('Warning! Error with midrise calculation!')
                            end
                        end
                    end
                end

                for ii = 1:S_E_L
                    WAVG{1,jj}{1,ci}(ii,kk) = nanmean(wavgs(ii,:)); %per experiment, per cell, per trace
                end
                
                if figures_on == 1
                    figure('position',[56 200 1400 490])
                    subplot(1,3,3)
                    hold on
                    title('waveform average')
                    for ii = 1:numel(timeindx)
                        if isfinite(nanmean(wavgs(:,ii)))==1
                            plot(wavgs(:,ii),'b')
                        end
                    end
                    plot(WAVG{1,jj}{1,ci}(:,kk),'r','LineWidth',3);
                    shg
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                mean_raw_amp = nanmean(RAW_AMP);
                total_time = numel(~isnan(aDAT))/samplert; % in second
                counts = sum(~isnan(sm_pk_ind)); 
                frequency = counts/total_time; % in Hz

                %%%% All results are saved as the following: 
                %%%% Per experiment (cell array)-> Per cell (cell array)-> Per
                %%%% trace (each row in the first column)
                FRQ{1,jj}{1,ci}(kk,1) = frequency;
                AMP_ALL{1,jj}{1,ci}(1:numel(RAW_AMP),kk) = RAW_AMP;

                meanAmp{1,jj}{1,ci}(kk,1) = mean_raw_amp;
                mini_count{1,jj}{1,ci}(kk,1) = counts;
                TOTAL_TIME{1,jj}{1,ci}(kk,1) = total_time;
                
                RISE{1,jj}{1,ci}(1:numel(RISE_TIME),kk) = RISE_TIME;
                DECAY{1,jj}{1,ci}(1:numel(DECAY_TIME),kk) = DECAY_TIME;
                wavgs_raw{1,jj}{1,ci}{kk} = wavgs;
                TIME_INDICES{1,jj}{1,ci}(1:numel(timeindx),kk) = timeindx;
                aDAT_all{1,jj}{1,ci}(1:numel(aDAT),kk) = aDAT;

                EVENT_START_INDICES{1,jj}{1,ci}(1:numel(event_start_ind),kk) = event_start_ind;
                PT3_AVERAGES{1,jj}{1,ci}{kk} = pt3_avg;
                SMOOTH_PEAK_INDICES{1,jj}{1,ci}(1:numel(sm_pk_ind),kk) = sm_pk_ind;
                DECAY_INDICES{1,jj}{1,ci}(1:numel(decay_ind),kk) = decay_ind;
                PT3_DECAY_VALUES{1,jj}{1,ci}{kk} = pt3_decay;

                %remove zeros from some cell arrays for plotting
                temp_1 = nonzeros(AMP_ALL{1,jj}{1,ci}(:,kk));
                AMP_ALL{1,jj}{1,ci}(:,kk) = NaN;
                AMP_ALL{1,jj}{1,ci}(1:numel(temp_1),kk) = temp_1;

                temp_2 = nonzeros(RISE{1,jj}{1,ci}(:,kk));
                RISE{1,jj}{1,ci}(:,kk) = NaN;
                RISE{1,jj}{1,ci}(1:numel(temp_2),kk) = temp_2;

                temp_3 = nonzeros(DECAY{1,jj}{1,ci}(:,kk));
                DECAY{1,jj}{1,ci}(:,kk) = NaN;
                DECAY{1,jj}{1,ci}(1:numel(temp_3),kk) = temp_3;

                temp_4 = nonzeros(TIME_INDICES{1,jj}{1,ci}(:,kk));
                TIME_INDICES{1,jj}{1,ci}(:,kk) = NaN;
                TIME_INDICES{1,jj}{1,ci}(1:numel(temp_4),kk) = temp_4;

                temp_5 = nonzeros(aDAT_all{1,jj}{1,ci}(:,kk));
                aDAT_all{1,jj}{1,ci}(:,kk) = NaN;
                aDAT_all{1,jj}{1,ci}(1:numel(temp_5),kk) = temp_5;

                temp_6 = nonzeros(SMOOTH_PEAK_INDICES{1,jj}{1,ci}(:,kk));
                SMOOTH_PEAK_INDICES{1,jj}{1,ci}(:,kk) = NaN;
                SMOOTH_PEAK_INDICES{1,jj}{1,ci}(1:numel(temp_6),kk) = temp_6;

                temp_7 = nonzeros(DECAY_INDICES{1,jj}{1,ci}(:,kk));
                DECAY_INDICES{1,jj}{1,ci}(:,kk) = NaN;
                DECAY_INDICES{1,jj}{1,ci}(1:numel(temp_7),kk) = temp_7;


                %%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                if figures_on == 1 
                    disp('plotting results:')

                    subplot(1,3,1:2)
                    hold on
                    title(strcat('Cell',num2str(ci),': trace',sname),'Interpreter', 'none')
                    plot(aDAT_all{1,jj}{1,ci}(:,kk),'k')


                    timeindx_slopeBL = timeindx_slopeBL(isnan(timeindx_slopeBL)==0);
                    timeindx_rise = timeindx_rise(isnan(timeindx_rise)==0);
                    timeindx_decay = timeindx_decay(isnan(timeindx_decay)==0);
                    timeindx_charge = timeindx_charge(isnan(timeindx_charge)==0);
                    timeindx_BLBAD = timeindx_BLBAD(isnan(timeindx_BLBAD)==0);
                    timeindx_noiseSTD = timeindx_noiseSTD(isnan(timeindx_noiseSTD)==0);
                    timeindx_multipk = timeindx_multipk(isnan(timeindx_multipk)==0);
                    timeindx_PeakAtEnd = timeindx_PeakAtEnd(isnan(timeindx_PeakAtEnd)==0);
                    timeindx_StartAfterPeak = timeindx_StartAfterPeak(isnan(timeindx_StartAfterPeak)==0);
                    timeindx_LatePeak = timeindx_LatePeak(isnan(timeindx_LatePeak)==0);
                    timeindx_decay_start = timeindx_decay_start(isnan(timeindx_decay_start)==0);
                    timeindx_refractory_per = timeindx_refractory_per(isnan(timeindx_refractory_per)==0);
                    timeindx_decayfit = timeindx_decayfit(isnan(timeindx_decayfit)==0);
                    timeindx_BLnopts =  timeindx_BLnopts(isnan(timeindx_BLnopts)==0);
                    timeindx_no_decay = timeindx_no_decay(isnan(timeindx_no_decay)==0);
                    timeindx_ampcutoff = timeindx_ampcutoff(isnan(timeindx_ampcutoff)==0);


                    scatter(timeindx,repmat( 15 ,1,length(timeindx)),'ro','filled'); % time index of each event
                    scatter(timeindx_slopeBL,repmat( 15.5,1,length(timeindx_slopeBL)),'bx') % baseline slope
                    scatter(timeindx_rise,repmat( 16,1,length(timeindx_rise)),'gx') % rise time
                    scatter(timeindx_decay,repmat( 16.5,1,length(timeindx_decay)),'rx') %decay shorter than rise
                    scatter(timeindx_charge,repmat( 17,1,length(timeindx_charge)),'bo','filled') %charge 
                    scatter(timeindx_BLBAD,repmat( 17.5,1,length(timeindx_BLBAD)),'ko','filled') % bad baseline
                    scatter(timeindx_noiseSTD,repmat( 18,1,length(timeindx_noiseSTD)),'mx') %std of noise
                    scatter(timeindx_multipk,repmat( 18.5,1,length(timeindx_multipk)),'cx') %mulitple peaks
                    scatter(timeindx_PeakAtEnd,repmat( 19,1,length(timeindx_PeakAtEnd)),'r*') %peak matches the last value in 3-pts average
                    scatter(timeindx_StartAfterPeak,repmat( 19.5,1,length(timeindx_StartAfterPeak)),'y*') %event start after peak
                    scatter(timeindx_LatePeak,repmat( 20.5,1,length(timeindx_LatePeak)),'b*') %event starts too late
                    scatter(timeindx_decay_start,repmat( 21,1,length(timeindx_decay_start)),'b+')%can't find decay start point
                    scatter(timeindx_refractory_per,repmat( 21.5,1,length(timeindx_refractory_per)),'g+') %time btw events too short
                    scatter(timeindx_decayfit,repmat(22,1,length(timeindx_decayfit)),'r+') %decay is not a good exp fit
                    scatter(timeindx_BLnopts,repmat(21.25,1,length(timeindx_BLnopts)),'k+') %not enought points for baseline
                    scatter(timeindx_no_decay,repmat(21.75,1,length(timeindx_no_decay)),'mo','filled') %decay cannot be found
                    scatter(timeindx_ampcutoff,repmat(20.75,1,length(timeindx_ampcutoff)),'co','filled') %amp off bounds
                    
                    for i = 1:numel(timeindx)
                        if isnan(RAW_AMP(i)) == 0
                            plot(event_start_ind{i}+timeindx(i)-1,nanmean(aDAT(event_start_ind{i}+...,
                                timeindx(i)-numptsBL-1:event_start_ind{i}+timeindx(i)-1)),'go') % event starts
                            plot(sm_pk_ind(i)+timeindx(i)-1,pt3_avg(sm_pk_ind(i),i),'bo') % peak position
                            if isnan(decay_ind(i)) == 0
                                plot(decay_ind(i)+timeindx(i)+sm_pk_ind(i),pt3_decay(decay_ind(i),i),'ro') % event ends
                            end
                        end
                    end
                end %%%%plot results

                s1 = 'Finished: ';
                disp(strcat(s1,'cell', num2str(ci), ' trace', num2str(kk)))

                fclose('all');
               end
            end %%%%file
        end
    end %%%%cell
end %%%% experiment

%%
%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = input('Save results? (1 = y, 2 = n):','s');
result = str2double(result);

if result == 1
    if risetime_cutoff >5 %(used 1.5 as acutal cutoff because too many excellent events would be excluded if set at 1)
        cd(strcat(fp_analyzed_data,'mIPSC'))
    end
    
        save(save_file_name,'AMP_ALL','WAVG','FRQ','wavgs_raw',...
            'RISE','TIME_INDICES',...
            'aDAT_all','EVENT_START_INDICES','PT3_AVERAGES','SMOOTH_PEAK_INDICES',...
            'DECAY_INDICES','PT3_DECAY_VALUES','TOTAL_TIME', 'HALF_WIDTH', 'DECAY',...
            'mini_count','meanAmp')
    
end
