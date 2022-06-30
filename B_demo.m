clc; clear; close all; 
addpath(genpath('matlab-src'));

%% Select demo file name here
data_file = 'data/GW_20150806_001.mat'; 
load(data_file, 'recordings');

%% Get data out 
ind_trace = 16; % demo sweep number 
len_trace = length(recordings(1).data); % number of data points 
dt = recordings(1).dt * 1e3; % seconds to ms, cuz whole-cell 
Vm = recordings(ind_trace).data * 1e3; % V to mV, cuz whole-cell 
t = (0:len_trace-1)*dt; % time vec 

%% Stimulus info 
stim_info = struct; 
stim_info.pulse_on_time         = 100; % time current injection on [ms] 
stim_info.pulse_off_time        = 600; % time current injection off [ms] 
stim_info.hyperpol_on_time      = 1000; % hyperpolarization on to determine Rin [ms]
stim_info.hyperpol_off_time     = 1100; % hyperpolarization off to determine Rin [ms]
stim_info.hyperpol_Iinj_amp     = 100; % amplitude of hyperpolarization current [pA]

%% Defining parameters 
spike_params = struct; 
spike_params.peak_prom          = 50; % peak prominence to detect AP, used in `findpeaks` [mV]
spike_params.min_dist           = 1; % minimal distance beteween peaks [ms]
spike_params.min_height         = 0; % minimum AP absolute heights, sometimes AP can go lower than 0mV [mV]
spike_params.min_width          = 0.05; % minumum width of AP, sometimes high amplitude but sharp noise is detected [ms]
spike_params.analyze_in_stim    = true; % just leave this on, too lazy to explain 

AP_params = struct; 
AP_params.max_tpre              = 1; % max time to look backward to find threshold [ms]
AP_params.min_dV_dt             = 8; % min dV/dt to detect AP thres [mV/ms or units]
AP_params.scale_min_dV_dt       = false; % if true, scale `min_dV_dt` by `max(dV/dt)`, to find threshold
% AP_params.min_dV_dt             = 0.1; % comment out if wanting to use the scaled version instead
% AP_params.scale_min_dV_dt       = true; 
AP_params.min_pts_cross_thres   = 2; % min # points crossing the set `min_dV_dt` to qualify as threshold
AP_params.max_tpost             = 2; % max time to look forward (could just use the next AP time) to find fAHP
AP_params.rise_start_APfactor   = 0.1; % min AP factor (between threshold and AP amplitude) to consider rise time
AP_params.rise_stop_APfactor    = 0.9; % max AP factor (between threshold and AP amplitude) to consider rise time
AP_params.rise_interp_factor    = 100; % factor to interp rise (between threshold and AP amplitude) 
AP_params.fall_start_APfactor   = 0.9; % max AP factor (between AP amplitude and minimum/fAHP amplitude) to consider fall time
AP_params.fall_stop_APfactor    = 0.1; %  min AP factor (between AP amplitude and minimum/fAHP amplitude) to consider fall time
AP_params.fall_interp_factor    = 100; % factor to interp rise (between AP amplitude and minimum/fAHP amplitude)

mAHP_params = struct; 

% find_opt: which option to use to estimate AHP parameters 
% - 'raw': get min from raw 
% - 'smooth': smoothed then min
% - 'doubexp': fitted with double exp then min
% - 'lowpass': low pass then min 
% - 'all': return all 
mAHP_params.find_opt            = 'all'; 
mAHP_params.max_tpost           = 300; % max time to look forward [ms]
mAHP_params.save_processed      = true; % whether to save the processed
mAHP_params.smooth_method       = 'movmean'; % to use when `find_opt='smooth'`
mAHP_params.smooth_tau          = 5; % to use when `find_opt='smooth'`
mAHP_params.lowpass_freq        = 50; % to use when `find_opt='lowpass'`

Vrest_calc_time_window          = 20; % window before stim on to calc rest
Rin_calc_time_window            = 20; % window before and after hyperpol to calc Rin

%% Analyze
ind_start_mAHP = ceil((stim_info.pulse_off_time)/dt);

spike_properties = calculate_spike_properties(t, Vm, spike_params, stim_info);
AP_properies = arrayfun(@(AP_ind) calculate_AP_properties(Vm, AP_ind, dt, AP_params), ...
    spike_properties.spike_inds, 'uni', 1);
AP_properies = structarray_to_struct(AP_properies);

mAHP_properties = calculate_mAHP(Vm, ind_start_mAHP, dt, mAHP_params);
mAHP_methods = fieldnames(mAHP_properties);

Vrest_properties = calculate_Vrest(t, Vm, stim_info, Vrest_calc_time_window);
Rin_properties = calculate_Rin(t, Vm, stim_info, Rin_calc_time_window);

%% Plot
figure; hold on; 

% 1. Plot membrane potential
plot(t, Vm, 'color', [0.8,0.8,0.8], 'displayname', 'raw', 'tag', 'legendon');

% 2. Plot AP locations
plot(t(spike_properties.spike_inds), Vm(spike_properties.spike_inds), 'o', 'displayname', 'AP', 'tag', 'legendon');

% 3. Plot spiking thresholds 
plot(AP_properies.AP_thres_time, AP_properies.AP_thres_Vm, 'o', 'displayname', 'thres', 'tag', 'legendon');

% 4. Plot fast AHP estimation 
plot(AP_properies.fAHP_time, AP_properies.fAHP_Vm, 'o', 'displayname', 'fAHP', 'tag', 'legendon');

% 5. Plot AP rise and fall time estimations
plot([AP_properies.AP_rise_start_time, AP_properies.AP_rise_stop_time]', ...
    [AP_properies.AP_rise_start_Vm,  AP_properies.AP_rise_stop_Vm]', '-r', ...
    'displayname', 'rise')
plot([AP_properies.AP_fall_start_time, AP_properies.AP_fall_stop_time]', ...
    [AP_properies.AP_fall_start_Vm,  AP_properies.AP_fall_stop_Vm]', '-b', ...
    'displayname', 'fall')

% 6. Plot AP widths 
plot([AP_properies.AP_width_start_time, AP_properies.AP_width_stop_time]', ...
    [AP_properies.AP_width_Vm,  AP_properies.AP_width_Vm]', '-k', ...
    'displayname', 'AP width')

% 7. Plot mAHP estimations from different methods
cmap = lines(length(mAHP_methods));
for i = 1:length(mAHP_methods)
    mAHP_meth = mAHP_methods{i}; 
    mAHP_prop = mAHP_properties.(mAHP_meth);
    if i > 1
    plot(mAHP_prop.processed_t, mAHP_prop.processed_Vm, ':', 'color', [cmap(i,:), 0.5], ...
        'displayname', ['find mAHP with ' mAHP_meth], 'tag', 'legendon'); 
    end
    plot(mAHP_prop.mAHP_time, mAHP_prop.mAHP_amp, 's', 'markersize', 10, ...
        'color', cmap(i,:), 'displayname', ['mAHP - ' mAHP_meth], 'tag', 'legendon'); 
end

xlabel('time (ms)');
ylabel('V_m (mV)');

legend(findobj(gca, 'tag', 'legendon'));

