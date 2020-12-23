load('20201201_baseline.mat', 'recordings');
%%
n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt; 
dt_ms = dt*1e3;
t = (0:len_trace-1)*dt_ms;

figure; hold on;
for i = 1:n_trace
plot(t, recordings(i).data*1e3 + 20*i, '-k')
end

%%
peak_reqs = {'MinPeakProminence', 50, 'MinPeakDistance', 1};
Vm = recordings(1).data*1e3; 
figure; hold on;
plot(t, Vm, '-k');

% findpeaks(Vm, t, 'MinPeakProminence', 50, 'MinPeakDistance', 1);
[~, pk_locs] = findpeaks(Vm, 'MinPeakProminence', 50, 'MinPeakDistance', 1/dt_ms);
plot(t(pk_locs), Vm(pk_locs), 'o');

max_tpre_ms = 1; 
thres_dV_dt = 5; % V/sec


ind_spkthres = arrayfun(@(x) ...
    firing_threshold(Vm, x, dt_ms, max_tpre_ms, thres_dV_dt), pk_locs);

plot(t(ind_spkthres), Vm(ind_spkthres), 'sb', 'markersize', 10);


yyaxis right;
plot(t(1:end-1), diff(Vm)/dt_ms);

%%
figure; hold on;
% plot(t, Vm, '-k');

findpeaks(-Vm, t, 'MinPeakProminence', 10, 'MinPeakDistance', 1);
% [~, pk_locs] = findpeaks(Vm, 'MinPeakProminence', 50, 'MinPeakDistance', 1/dt_ms);
% plot(t(pk_locs), Vm(pk_locs), 'o');
%%

% for detecting spikes (APs)
spk_peak_prom       = 50; % mV, spike prominence
spk_min_dist        = 1; % ms, min spike distances
spk_min_height      = 0; % mV, min spike amplitude 

% for detecting spike thresholds
max_tpre_ms         = 1; % ms, max time to go backwards relative to spike
thres_dV_dt         = 5; % V/sec, threshold to detect "kink" point
min_pts_cross_thres = 2; % min # of points that have to cross `thres_dV_dt`

figure; hold on; 
for i = 1:n_trace
    Vm = recordings(i).data * 1e3; 
    
    [~, pk_locs] = findpeaks(Vm, 'MinPeakProminence', spk_peak_prom, ...
        'MinPeakDistance', spk_min_dist/dt_ms, 'MinPeakHeight', spk_min_height);
    
    ind_spkthres = arrayfun(@(x) firing_threshold(Vm, x, dt_ms, ...
        max_tpre_ms, thres_dV_dt, min_pts_cross_thres), pk_locs);
    
    plot(t, Vm, '-k');
    plot(t(pk_locs), Vm(pk_locs), 'o');
    plot(t(ind_spkthres), Vm(ind_spkthres), 'sb', 'markersize', 10);
    title(sprintf('trace %02d', i)); 
    
    xlim([100,700]); %ylim([-60, -20]);
    waitforbuttonpress
    cla;
end

