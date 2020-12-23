clc; clear; close all; 
run matlab_startup.m; 

%% Change file name here
% file_name = 'gw150806_2_Baseline'; 
file_name = 'gw150806_2_Post';
data_path = 'sample-data/20150806'; 

% file_name = '20200724_1_baseline'; 
% data_path = 'sample-data/20200724'; 

load(fullfile(data_path, [file_name '.mat']), 'recordings');

%% Get data out 
n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt * 1e3;
Vm = recordings(106-31).data * 1e3; 
t = (0:len_trace-1)*dt;

%% Stimulus info 
stim_info = struct; 
stim_info.pulse_on_time         = 100; % ms 
stim_info.pulse_off_time        = 600; % ms 
stim_info.hyperpol_on_time      = 1000; % ms
stim_info.hyperpol_off_time     = 1100; % ms
stim_info.hyperpol_Iinj_amp     = 100; % pA

%% Defining parameters 
spike_params = struct; 
spike_params.peak_prom          = 50; 
spike_params.min_dist           = 1; 
spike_params.min_height         = 0; 
spike_params.min_width          = 0.05;  
spike_params.analyze_in_stim    = true; 

AP_params = struct; 
AP_params.max_tpre              = 1; 
AP_params.min_dV_dt             = 8;  
AP_params.scale_min_dV_dt       = false; 
% AP_params.min_dV_dt             = 0.1;  
% AP_params.scale_min_dV_dt       = true; 
AP_params.min_pts_cross_thres   = 2; 
AP_params.max_tpost             = 0.5; 
AP_params.rise_start_APfactor   = 0.1;
AP_params.rise_stop_APfactor    = 0.9;
AP_params.rise_interp_factor    = 100;
AP_params.fall_start_APfactor   = 0.9;
AP_params.fall_stop_APfactor    = 0.1; 
AP_params.fall_interp_factor    = 100;

mAHP_params = struct; 
mAHP_params.find_opt            = 'all'; 
mAHP_params.max_tpost           = 300; 
mAHP_params.save_processed      = true; 
mAHP_params.smooth_method       = 'movmean'; 
mAHP_params.smooth_tau          = 5; 
% mAHP_params.doubexp_lowbound    = [-50,-100,0,0,-90]; 
% mAHP_params.doubexp_uppbound    = [100,50,150,500,-40];
% mAHP_params.doubexp_startpnt    = [10,-10,50,100,-60];
mAHP_params.lowpass_freq        = 50; 

Vrest_calc_time_window          = 20; 
RIn_calc_time_window            = 20; 

%% Analyze
ind_start_mAHP = ceil((stim_info.pulse_off_time)/dt);

spike_properties = calculate_spike_properties(t, Vm, spike_params, stim_info);
AP_properies = arrayfun(@(AP_ind) calculate_AP_properties(Vm, AP_ind, dt, AP_params), ...
    spike_properties.spike_inds, 'uni', 1);
AP_properies = structarray_to_struct(AP_properies);

mAHP_properties = calculate_mAHP(Vm, ind_start_mAHP, dt, mAHP_params);
mAHP_methods = fieldnames(mAHP_properties);

Vrest_properties = calculate_Vrest(t, Vm, stim_info, Vrest_calc_time_window);
Rin_properties = calculate_Rin(t, Vm, stim_info, RIn_calc_time_window);

%% Plot
figure; hold on; 
plot(t, Vm, 'color', [0.8,0.8,0.8], 'displayname', 'raw');
plot(t(spike_properties.spike_inds), Vm(spike_properties.spike_inds), 'o', 'displayname', 'AP');

plot(AP_properies.AP_thres_time, AP_properies.AP_thres_Vm, 'o', 'displayname', 'thres');
plot(AP_properies.fAHP_time, AP_properies.fAHP_Vm, 'o', 'displayname', 'fAHP');

plot([AP_properies.AP_rise_start_time, AP_properies.AP_rise_stop_time]', ...
    [AP_properies.AP_rise_start_Vm,  AP_properies.AP_rise_stop_Vm]', '-r', ...
    'displayname', 'rise')

plot([AP_properies.AP_fall_start_time, AP_properies.AP_fall_stop_time]', ...
    [AP_properies.AP_fall_start_Vm,  AP_properies.AP_fall_stop_Vm]', '-b', ...
    'displayname', 'fall')

plot([AP_properies.AP_width_start_time, AP_properies.AP_width_stop_time]', ...
    [AP_properies.AP_width_Vm,  AP_properies.AP_width_Vm]', '-k', ...
    'displayname', 'width')

cmap = prism(length(mAHP_methods));
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


legend(findobj(gca, 'tag', 'legendon'));


%%
file_name = 'gw150806_2_Post'; 
data_path = 'sample-data/20150806';
load(fullfile(data_path, [file_name '.mat']), 'recordings');


n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt * 1e3;
t = (0:len_trace-1)*dt;

saved_measures = cell(n_trace,1); 

for i = 1:n_trace
    Vm = recordings(i).data * 1e3; 
    
    spike_properties = calculate_spike_properties(t, Vm, spike_params, stim_info);
    if spike_properties.num_spikes > 0 
    AP_properies = arrayfun(@(AP_ind) calculate_AP_properties(Vm, AP_ind, dt, AP_params), ...
        spike_properties.spike_inds, 'uni', 1);
    
    AP_properies = structfun(@(x) nanmean(x), structarray_to_struct(AP_properies), 'uni', 0);
    else
        fn_APs = fieldnames(AP_properies);
        for j = 1:length(fn_APs)
            AP_properies.(fn_APs{j}) = nan;
        end
    end
    mAHP_properties = calculate_mAHP(Vm, ind_start_mAHP, dt, mAHP_params);
    mAHP_properties = rmfield(mAHP_properties, 'raw'); 
    mAHP_amp = mean(structfun(@(x) x.mAHP_amp, mAHP_properties));
    mAHP_time = mean(structfun(@(x) x.mAHP_time, mAHP_properties));
    
    Vrest_properties = calculate_Vrest(t, Vm, stim_info, Vrest_calc_time_window);
    Rin_properties = calculate_Rin(t, Vm, stim_info, RIn_calc_time_window);
    
    tmp_struct = struct; 
    tmp_struct.Rin = Rin_properties.Rin; 
    tmp_struct.Vrest = Vrest_properties.Vrest; 
    tmp_struct.mAHP_amp = mAHP_amp;
    tmp_struct.mAHP_time = mAHP_time;
    tmp_struct.firing_rate = spike_properties.firing_rate;
    tmp_struct.latency = spike_properties.latency;
    tmp_struct.CV_ISI = spike_properties.CV_ISI;
    tmp_struct.mean_ISI = spike_properties.mean_ISI;
    tmp_struct.adapt_index = spike_properties.adapt_index;
    
    saved_measures{i} = mergefield_struct(tmp_struct, AP_properies); 
    i
end
saved_measures = structarray_to_struct(horzcat(saved_measures{:})); 
%%
figure; 
plot(saved_measures.Vrest)

% base_measures = saved_measures;
% post_measures = saved_measures;
%%
all_fields = fieldnames(base_measures); 
all_measures = struct;
for i = 1:length(all_fields)
    f_i = all_fields{i};
    all_measures.(f_i) = horzcat(base_measures.(f_i), post_measures.(f_i));
end
%%
figure; 
plot(all_measures.AP_thres_Vm - all_measures.Vrest);
