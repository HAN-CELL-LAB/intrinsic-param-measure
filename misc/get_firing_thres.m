clc; clear; close all; 

%% Change file name here
file_name = '20201201_baseline'; 

load([file_name '.mat'], 'recordings');

%% Get data out 
n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt; 
dt_ms = dt*1e3;
t = (0:len_trace-1)*dt_ms;

%% Set parameters for finding spikes and thresholds
% change until satisfied

% for detecting spikes (APs)
spk_peak_prom       = 50; % mV, spike prominence (30-50 mV is usually good)
spk_min_dist        = 1; % ms, min spike distances (1-3 ms is good)
spk_min_height      = 0; % mV, min spike amplitude (0mV would be good, unless there are spikes under 0mV, then maybe -10mV)
spk_min_width       = 0.05; % ms, min spike width in order to filter out abrupt artifact 

% for detecting spike thresholds
max_tpre_ms         = 1; % ms, max time to go backwards relative to spike (0.5 to 1 or 2 ms would be good)
thres_dV_dt         = 5; % V/sec, threshold to detect "kink" point (5, 10 or 20 would be good) 
min_pts_cross_thres = 2; % min # of points that have to cross `thres_dV_dt` (2-3 would be good)

%% Change above and run this section to inspect if everything's good
% if not, change parameters above 

spike_threshold_data = cell(n_trace,1); 
spike_threshold_time = cell(n_trace,1); 
spike_threshold_inds = cell(n_trace,1); 

figure; hold on; 
for i = 1:n_trace
    Vm = recordings(i).data * 1e3; 
    
    [~, pk_locs] = findpeaks(Vm, 'MinPeakProminence', spk_peak_prom, ...
        'MinPeakDistance', spk_min_dist/dt_ms, 'MinPeakHeight', spk_min_height, ...
        'MinPeakWidth', spk_min_width/dt_ms);
    
    ind_spkthres = arrayfun(@(x) firing_threshold(Vm, x, dt_ms, ...
        max_tpre_ms, thres_dV_dt, min_pts_cross_thres), pk_locs);
    
    AP_ind = pk_locs;
    AP_Vm = Vm(pk_locs);
    AP_t = t(pk_locs);
    AP_n = length(pk_locs); 
    
    
    thres_ind = ind_spkthres;
    thres_t = t(ind_spkthres);
    thres_Vm = Vm(ind_spkthres);
    
    plot(t, Vm, '-k');
    plot(AP_t, AP_Vm, 'o');
    plot(thres_t, thres_Vm, 'sb', 'markersize', 10);
    title(sprintf('trace %02d', i)); 
    
    spike_threshold_inds{i} = ind_spkthres;
    spike_threshold_time{i} = t(ind_spkthres);
    spike_threshold_data{i} = Vm(ind_spkthres); 
    
    xlim([100,700]); %ylim([-60, -20]);
    waitforbuttonpress
    cla;
end

%% If above is ok, run this section to save 
params = struct;
params.spk_peak_prom = spk_peak_prom;
params.spk_min_dist = spk_min_dist;
params.spk_min_height = spk_min_height;
params.max_tpre_ms = max_tpre_ms;
params.thres_dV_dt = thres_dV_dt;
params.min_pts_cross_thres = min_pts_cross_thres;

thres_file_name = [file_name '-firingthresholds']; 
save([thres_file_name '.mat'], 'params', 'spike_threshold_data', ...
    'spike_threshold_inds', 'spike_threshold_time'); 

writetable(cell2table(spike_threshold_data), [thres_file_name '.csv'])

%%
    
fit_param_list = cell(n_trace,1);
fit_perf_list = cell(n_trace,1);
figure; hold on; 
for i = 1:n_trace
    Vm = recordings(i).data * 1e3; 
    

t_e = t(t>602 & t <900) - 602;
V_e = Vm(t>602 & t <900);

fit_opt_obj = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10,-80,0,0,-90],...
               'Upper',[80,10,100,300,-40],...
               'StartPoint',[10,-10,10,50,-60]);
fit_type_obj = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2) + y_c', 'options', fit_opt_obj);
[fit_param, fit_perf] = fit(t_e',V_e,fit_type_obj);
ft_V_e = fit_param(t_e);

fit_param_list{i} = fit_param; 
fit_perf_list{i} = fit_perf;
tau_smoothwin = 5; 
sm_V_e = smoothdata(V_e, 'movmean', ceil(tau_smoothwin/dt_ms)); 

[~, ind_min_raw] = min(V_e); 
[~, ind_min_fit] = min(ft_V_e); 
[~, ind_min_smooth] = min(sm_V_e); 

plot(t_e, V_e, '-', 'color', [0.7,0.7,0.7], 'displayname', 'raw'); 
plot(t_e(ind_min_raw), V_e(ind_min_raw), 'o', 'color', [0.7,0.7,0.7], 'displayname', 'raw min', 'markersize', 20);
plot(t_e, ft_V_e, ':r', 'displayname', 'fit double-exp');
plot(t_e(ind_min_fit), ft_V_e(ind_min_fit), 'ro', 'displayname', 'min of fit', 'markersize', 20);

tau_smoothwins = [2,10]; 
c_taus = [0.4,0.4,0.8; 0,0,1]; 

   title(sprintf('trace %02d', i)); 
   
for j = 1:length(tau_smoothwins)
    tau_smoothwin = tau_smoothwins(j);
    sm_V_e = smoothdata(V_e, 'movmean', ceil(tau_smoothwin/dt_ms));
    [~, ind_min_smooth] = min(sm_V_e);
    
    plot(t_e, sm_V_e, '-', 'color',  c_taus(j,:), 'displayname', sprintf('smoothed %gms', tau_smoothwin));
    plot(t_e(ind_min_smooth), sm_V_e(ind_min_smooth), 'o', 'color',  c_taus(j,:), ...
        'displayname', sprintf('min of smooth %gms', tau_smoothwin), 'markersize', 20);
end

legend('show');

    waitforbuttonpress
    cla;
end
%%
fit_param_list = cellfun(@(x) struct('a1', x.a1, 'a2', x.a2, ...
    'tau1', x.tau1, 'tau2', x.tau2, 'y_c', x.y_c), fit_param_list, 'uni', 1);
fit_perf_list = vertcat(fit_perf_list{:});
%%
% mahp_tau1s = cellfun(@(x) x.tau1, 
figure;
subplot(311); hold on;
plot([fit_param_list.a1], 'o');
plot([fit_param_list.a2], 'o');
ylabel('$a_1$ or $a_2$');
yyaxis right; 
plot([fit_param_list.y_c], 'o');
ylabel('$y_c$');
legend('$a_1$','$a_2$', '$y_c$');
title('$AHP = a_1*\exp(-t/\tau_1) + a_2*\exp(-t/\tau_2) + y_c$');

subplot(312); hold on;
plot([fit_param_list.tau1], 'o');
plot([fit_param_list.tau2], 'o');
legend('$\tau_1$','$\tau_2$');
ylabel('$\tau_1$ or $\tau_2$');

subplot(313); hold on;
plot([fit_perf_list.rsquare], 'o');
ylabel('$R^2$');
yyaxis right
plot([fit_perf_list.rmse], 'o');
ylabel('RMSE');
legend('$R^2$', 'RMSE');

%%
figure; 
subplot(311); hold on; 
histogram([fit_param_list.a1],50, 'BinLimits', [-100,100])
histogram([fit_param_list.a2],50, 'BinLimits', [-100,100])
histogram([fit_param_list.y_c],50, 'BinLimits', [-100,-50])
legend('$a_1$','$a_2$','$y_c$');

subplot(312); hold on;
histogram([fit_param_list.tau1],100, 'BinLimits', [0, 500])
histogram([fit_param_list.tau2],100, 'BinLimits', [0, 500])
legend('$\tau_1$','$\tau_2$');

subplot(313); hold on;
histogram([fit_perf_list.rsquare],40);
legend('$R^2$');