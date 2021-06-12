%% import mat files, areas of each variable stored under each experimental condition

%location where the mat file will be saved
fp_analyzed_data = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\area_under_curve_fI';

%name of the saved file
filename = 'Area_Under_Curve_wt_cp_L5_fI.mat';

%save results
save_file = 1;

% names of experimental conditions within the analyzed fI mat file
condn = {'WT_CNO_R','DR_CNO_R','WT_CNO_B','DR_CNO_B'};

%field names
fn = {'MFR', 'IFR', 'mean_IFR'}; 

%names of each group saved in the area mat file
%area_names = cell(1,numel(condn));
    
%% 
data_temp = cell(1,numel(condn));
dst = cell(1,numel(condn));
cell_num = NaN(numel(condn),1);

for gi = 1:numel(condn) %per condition
    %area_names{1,gi} = strcat(condn{1,gi},'_area');
    
   for fii = 1:numel(fn) %per field
       cell_num(gi,1) = size(eval(strcat(condn{gi},'.',fn{fii})),2);
       current_data = eval(strcat(condn{gi},'.',fn{fii}));
       
       for cii = 1:max(cell_num(gi)) %per cell
      
           if isnan(current_data(1:20,cii))
               data_temp{1,gi}{1,fii}(cii,1) = NaN;
           else
              data_temp{1,gi}{1,fii}(cii,1) = areaundercurve(current_inj,current_data(1:20,cii));
           end
          
       end
       
       dst{gi}.(fn{fii}) = data_temp{1,gi}{1,fii}(1:cell_num(gi,1));
       
   end
end

%save data to corresponding data structure (by experimental condition)
WT_CNO_R_area = dst{1};
DR_CNO_R_area = dst{2};
WT_CNO_B_area = dst{3};
DR_CNO_B_area = dst{4};
%% save file
if save_file == 1
    cd(fp_analyzed_data)
    
    condn_area = cell(1,numel(condn));
    for cdi = 1:numel(condn)
        condn_area{cdi} = strcat(condn{cdi},'_area');
    end
        
    save(filename,condn_area{1},condn_area{2},condn_area{3},condn_area{4})
end
