function res = calculate_spike_properties(t, Vm, params, stim_info, dry_run)

%% Initialization and exit if dry run
res = struct;

res.spike_inds = nan;
res.num_spikes = nan;
res.spike_times = nan;
res.firing_rate = nan;
res.latency = nan;
res.all_ISI = nan;
res.first_ISI = nan;
res.mean_ISI = nan;
res.std_ISI = nan;
res.CV_ISI = nan;
res.adapt_index = nan;

if ~exist('dry_run', 'var')
    dry_run = false;
end

if dry_run
    return;
end

%% Process parameters 
analyze_in_stim = get_field_or_default(params, 'analyze_in_stim', true); 
peak_prom       = get_field_or_default(params, 'peak_prom', 50);
min_dist        = get_field_or_default(params, 'min_dist', 1);
min_height      = get_field_or_default(params, 'min_height', 0);
min_width       = get_field_or_default(params, 'min_width', 0.05);

dt = t(2) - t(1); 
shift_ind = 0; 

if analyze_in_stim 
    if ~exist('stim_info', 'var') 
        error('When "analyze_in_stim" is true, need to provide `stim_info`');
    end
    pulse_on = stim_info.pulse_on_time;
    pulse_off = stim_info.pulse_off_time;
    shift_ind = find_nearest_ind(t, pulse_on) - 1;
    limit_time = t >= pulse_on & t <= pulse_off; 
    t = t(limit_time); 
    Vm = Vm(limit_time); 
end

[~, res.spike_inds] = findpeaks(Vm, ...
    'MinPeakProminence', peak_prom, 'MinPeakHeight', min_height, ...
    'MinPeakDistance', min_dist/dt, 'MinPeakWidth', min_width/dt);

res.num_spikes  = length(res.spike_inds); 
res.spike_times = t(res.spike_inds);
res.spike_inds  = res.spike_inds + shift_ind;

T_sec = (t(end)-t(1)) * 1e-3;
res.firing_rate = res.num_spikes/T_sec;

if res.num_spikes >= 1
    res.latency     = res.spike_times(1) - t(1);
else 
    res.latency     = nan;
end

res.all_ISI = diff(res.spike_times); 

if length(res.all_ISI) >= 1
    res.first_ISI = res.all_ISI(1);
    res.mean_ISI = mean(res.all_ISI);
else
    res.first_ISI = nan;
    res.mean_ISI = nan;
end

if length(res.all_ISI) >= 2
    res.std_ISI = std(res.all_ISI);
    res.CV_ISI = res.std_ISI / res.mean_ISI;
    res.adapt_index = mean(diff(res.all_ISI) ./ (res.all_ISI(1:end-1) + res.all_ISI(2:end)));
    
else
    res.std_ISI = nan;
    res.CV_ISI = nan;
    res.adapt_index = nan;
end

end