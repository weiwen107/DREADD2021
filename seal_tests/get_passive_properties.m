function [Rs_est_scaled, Rs_scaled, Rt_scaled, Cm_scaled, Vm_scaled] =...
    get_passive_properties(data, step_start, vstep, samprate, figure_on)

% data - a vector containing the trace with seal-test 
% step start - this is where the sealtest starts, should be in seconds
% vstep - the voltage of the step, in V
% the sampling rate, in Hz
% figures_on = 1 for plotting

data = data * 10^-12; % pA to A

V_test = vstep; 
V_hold = -7e-2; % holding potential (mV to V)

t=step_start+0.5; %starting point for fit, since a negative voltage is applied, starting point should be the end of the step.
%samprate = 10000; %in Hz
tracevals = data((t-0.2)*samprate:(t*samprate+6000)); %the range of data points to be fitted (in one voltage step)

peak = max(tracevals); 
startfit = find(tracevals==peak,1,'last');% offset to avoid the effects of pipette capacitance 

ss_Start = 400; %timepoints after the peak, used to define the steady state region
ss_End = 900; %timepoints after the peak, used to define the steady state region

pulse_end = 0.2*samprate;

baseline_val = tracevals(pulse_end - 1100:pulse_end - 100);
I_baseline_std = nanstd(baseline_val);
% for bi = 1:numel(baseline_val)
%     if baseline_val(bi)<(nanmean(baseline_val)-4*I_baseline_std)
%         baseline_val(bi) = NaN;
%     end
% end
I_baseline = nanmedian(baseline_val); %calculate the median baseline current (during the V step)

ss_val = tracevals(pulse_end + ss_Start:pulse_end + ss_End);
I_ss_std = nanstd(ss_val);
% for si = 1:numel(ss_val)
%     if ss_val(si)< nanmean(ss_val)-4*I_ss_std
%         ss_val(bi) = NaN;
%     end
% end
I_ss = nanmedian(ss_val); %calculate the median steady state current during the test

dI = I_ss - I_baseline; 

Rt = -V_test/dI; %input resistance 


Vm = -(I_ss * Rt) + V_hold; %membrane potential


whole_seal_fit_x = (startfit:startfit + ss_End)'; %from peak to the end of steady state
ss_seal_fit_x = (startfit + ss_Start:startfit + ss_End)';  
baseline = nanmedian(tracevals(ss_seal_fit_x)); %normalize the transient for fitting


transient_vals = tracevals(whole_seal_fit_x) - baseline; %normalized transient
transient_vals = transient_vals*10^12; %A to pA


transient_seal_x = (pulse_end+1:pulse_end + ss_Start)'; %both the rising and decay phase, without the steady state
transient_seal_vals = tracevals(transient_seal_x) - baseline;
transient_end_est = find(transient_seal_vals >= baseline, 1, 'last'); %estimated range of the area under the transient
Q1 = sum((abs(transient_seal_vals(1:transient_end_est).*(1/samprate)))); %sum of charge in the transient above the ss response 

tau_est = (find(transient_vals<(transient_vals(1)*0.37), 1,'first') - 1)/(samprate); %transient lifetime

%%%%%the fitting range of the transient
transient_st = find(transient_vals<(transient_vals(1)*0.95),1);
transient_end = find(transient_vals<(transient_vals(1)*0.05),1);

transient_vals = transient_vals(transient_st:transient_end);
%size(transient_vals)

time = (0:(1/samprate):(length(transient_vals)/samprate)-(1/samprate)).';

transient_vals = transient_vals - mean(transient_vals(end-5:end));
s = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [transient_vals(1), tau_est*0.5, transient_vals(7), tau_est*1.5],...
     'Lower', [mean(transient_vals(1:2)*.9), tau_est*0.001, mean(transient_vals(1:5))*.1, tau_est*0.01],...
     'Upper', [transient_vals(1)*1.1, tau_est*5, transient_vals(1)*1.25, tau_est*15]);
f = fittype('a*exp(-x/b) + c*exp(-x/d)','options',s);
%watchfithappen(time, transient_vals, f, 10)

[exp_fit,~] = fit(time,transient_vals,f);
cval = coeffvalues(exp_fit);

fast_exp = cval(1)*exp(-time./cval(2)); %values of the fast exponential component
slow_exp = cval(3)*exp(-time./cval(4)); %values of the slow exponential component
 
if figure_on == 1
     figure('position',[56 200 1200 490]); 
     subplot(1,4,1)
     plot(time, fast_exp); 
     hold on
     plot(time, slow_exp);   
end

%Q1_minus_slow = Q1 - sum(slow_exp.*(1/samprate))*10^-12; %%charge calculated by excluding the slow component

Q1_fromfit = sum(fast_exp.*(1/samprate))*10^-12; %charge calculated from first-exponential fit, usually lower than the estimated Q1


exp_Q2 = dI*exp(-time./cval(2));
Q2_fromfit = sum(exp_Q2.*(1/samprate));


% fast_Tau = cval(2);
% slow_Tau = cval(4);

%Q2 = dI * fast_Tau; %charges of the correction factor (of Q1)
Qt_calc = Q1_fromfit + Q2_fromfit; % total charges under the transient
Qt = Q1+Q2_fromfit;
%Qt_corr = Q1_minus_slow + Q2;
Cm = -Qt_calc/V_test; %membrane capacitance 
Cm_corr = -Qt/V_test;


num_pts_back_for_extrap = ((11+transient_st)-find(tracevals(startfit-10:startfit)>I_baseline+5*I_baseline_std,1));
%%%%need to compensate above line for the shift in the first transient
%%%%value above

%disp(num_pts_back_for_extrap)
[vals] = find(tracevals==peak);
if length(vals) < transient_st + 1 
  if num_pts_back_for_extrap > transient_st + 1 
      num_pts_back_for_extrap = transient_st + 1;
  end
end

Rs_est = -V_test/((exp_fit((-(num_pts_back_for_extrap/samprate)))*10^-12+I_baseline) + abs(I_baseline));

if isempty(Rs_est)
    Rs_est = NaN;
end

Rs = -V_test/(peak+abs(I_baseline));


if figure_on == 1
    subplot(1,4,2)
    plot(tracevals)
    subplot(1,4,3)
    plot(exp_fit,time,transient_vals)
    subplot(1,4,4)
    scatter(startfit-num_pts_back_for_extrap : startfit+transient_end, ...
        tracevals(startfit-num_pts_back_for_extrap : startfit+transient_end))
    hold on
    scatter(-num_pts_back_for_extrap+startfit,...
        ((exp_fit((-(num_pts_back_for_extrap/samprate)))*10^-12+I_baseline) + abs(I_baseline)),'r')
    %text(10,50,num2str(Rs*1e-6))
    %axis([0 50 min(transient_seal_vals) 100])
end

Rs_scaled = Rs*1e-6; %MOhm
Rs_est_scaled = Rs_est*1e-6; %MOhm
Rt_scaled = Rt*1e-6; %MOhm
Cm_scaled = Cm*1e12; %pF
Vm_scaled = Vm*1e3; %mV