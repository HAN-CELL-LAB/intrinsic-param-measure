function res = calculate_AP_properties(Vm, ind_AP, dt, params, dry_run)
% Obtain action potential (AP) properties, including AP amplitude,
% threshold, rise time, half-height width, fast AHP (minimum Vm after)
% INPUT:
%   Vm [mV]: vector of membrane potential
%   ind_AP [no unit]
%   dt [ms]
%   params: struct with
%       + max_tpre [ms]: max time to look backward to find threshold
%       + min_dV_dt [V/sec or mV/ms or units]: min dV/dt to detect AP thres 
%       + scale_min_dV_dt [boolean] [def=false]: if true, scale
%       `min_dV_dt` by `max(dV/dt)`, to find threshold
%       + min_pts_cross_thres [no unit] [def=3]: min # points crossing the
%       set `min_dV_dt` to qualify as threshold
%       + max_tpost [ms] [def=0.5] : max time to look forward (could just 
%       use the next AP time) to find fAHP
%       + rise_interp_factor [no unit] [def=100]: factor to interp rise
%       (between threshold and AP amplitude) 
%       + rise_start_APfactor [no unit] [def=0.1]: min AP factor (between
%       threshold and AP amplitude) to consider rise time
%       + rise_stop_APfactor [no unit] [def=0.9]: max AP factor (between
%       threshold and AP amplitude) to consider rise time
%           NOTE: default is to look for rise time by finding the time
%           from the 10% and the 90% between the threshold and AP amplitude
%       + fall_interp_factor [no unit] [def=100]: factor to interp rise
%       (between AP amplitude and minimum/fAHP amplitude)
%       + fall_start_APfactor [no unit] [def=0.9]: max AP factor (between
%       AP amplitude and minimum/fAHP amplitude) to consider fall time
%       + fall_stop_APfactor [no unit] [def=0.1]: min AP factor (between
%       AP amplitude and minimum/fAHP amplitude) to consider fall time
% OUTPUT: 
%   res: struct with 
% TODO: 
%   (1) edge case when spike happens at the start and end of trace (is it
%   necessary?)
%   (2) linear approx rise and fall time in addition to finding time via
%   indices nearest to compared Vm values 


%% Initialization and exit if dry run 
res = struct;

res.AP_Vm = nan;
res.AP_ind = nan;
res.AP_time = nan;
res.fAHP_Vm = nan;
res.fAHP_ind = nan;
res.fAHP_time = nan;
res.AP_thres_ind = nan;
res.AP_thres_time = nan;
res.AP_thres_Vm = nan;
res.AP_rise_start_time = nan;
res.AP_rise_stop_time = nan;
res.AP_rise_start_Vm = nan;
res.AP_rise_stop_Vm = nan;
res.AP_rise_time = nan;
res.AP_fall_start_time = nan;
res.AP_fall_stop_time = nan;
res.AP_fall_start_Vm = nan;
res.AP_fall_stop_Vm = nan;
res.AP_fall_time = nan;
res.AP_width_start_time = nan;
res.AP_width_stop_time = nan;
res.AP_width_Vm = nan;
res.AP_width = nan;

if ~exist('dry_run', 'var')
    dry_run = false;
end

if dry_run
    return; 
end

%% Processing parameters 
max_tpre    = params.max_tpre; 
min_dV_dt   = params.min_dV_dt; 

min_pts_cross_thres = get_field_or_default(params, 'min_pts_cross_thres', 3); 
scale_min_dV_dt = get_field_or_default(params, 'scale_min_dV_dt', false); 
max_tpost = get_field_or_default(params, 'max_tpost', 0.5); 


rise_interp_factor = get_field_or_default(params, 'rise_interp_factor', 100);
rise_start_APfactor = get_field_or_default(params, 'rise_start_APfactor', 0.1);
rise_stop_APfactor = get_field_or_default(params, 'rise_stop_APfactor', 0.9);

fall_interp_factor = get_field_or_default(params, 'fall_interp_factor', 100);
fall_start_APfactor = get_field_or_default(params, 'fall_start_APfactor', 0.9);
fall_stop_APfactor = get_field_or_default(params, 'fall_stop_APfactor', 0.1);

%% Find amplitude of AP
res.AP_Vm   = Vm(ind_AP);
res.AP_ind  = ind_AP;
res.AP_time = dt * (ind_AP - 1);

%% Find minimum (fast AHP)
incr_post = ceil(max_tpost/dt); 
[res.fAHP_Vm, fAHP_ind] = min(Vm(ind_AP:ind_AP+incr_post)); 
res.fAHP_ind = fAHP_ind + ind_AP - 1;
res.fAHP_time = dt * (res.fAHP_ind - 1); 

%% Finding threshold of AP
decr_pre = ceil(max_tpre/dt); 
Vpre = Vm(ind_AP-decr_pre:ind_AP);
dVpre_dt = diff(Vpre) / dt;  

if scale_min_dV_dt
    min_dV_dt = min_dV_dt * max(dVpre_dt); 
end

cmp_thres = dVpre_dt >= min_dV_dt;
ind_thres = find_start_of_sustained_on(cmp_thres, min_pts_cross_thres);

res.AP_thres_ind = return_nan_if_empty(ind_AP - decr_pre + ind_thres - 1); 
res.AP_thres_time = dt * (res.AP_thres_ind - 1); 

if ~isnan(res.AP_thres_ind)
    res.AP_thres_Vm  = Vm(res.AP_thres_ind);
else
    res.AP_thres_Vm = nan;
end



%% Find rise time of AP 
if ~isnan(res.AP_thres_ind)
    Vm_preAP = Vm(res.AP_thres_ind:res.AP_ind);
    thres_to_amp = Vm_preAP(end) - Vm_preAP(1);
    t_preAP = dt * (0:length(Vm_preAP)-1);
    t_preAP_interp = linspace(t_preAP(1), t_preAP(end), length(t_preAP)*rise_interp_factor); 
    
    Vm_preAP = interp1(t_preAP, Vm_preAP, t_preAP_interp); 
    shift_ind_pre = res.AP_thres_ind - 1;
    
    rise_start_Vm           = Vm_preAP(1) + rise_start_APfactor * thres_to_amp;
    rise_stop_Vm            = Vm_preAP(1) + rise_stop_APfactor * thres_to_amp;
    
    rise_start_ind          = find_nearest_ind(Vm_preAP, rise_start_Vm);
    rise_stop_ind           = find_nearest_ind(Vm_preAP, rise_stop_Vm);
    
    res.AP_rise_start_time  = t_preAP_interp(rise_start_ind) + shift_ind_pre*dt;
    res.AP_rise_stop_time   = t_preAP_interp(rise_stop_ind) + shift_ind_pre*dt;
    res.AP_rise_start_Vm    = Vm_preAP(rise_start_ind);
    res.AP_rise_stop_Vm     = Vm_preAP(rise_stop_ind);
    
    res.AP_rise_time        = abs(res.AP_rise_stop_time - res.AP_rise_start_time);
else
    res.AP_rise_start_time  = nan;
    res.AP_rise_stop_time   = nan;
    res.AP_rise_start_Vm    = nan;
    res.AP_rise_stop_Vm     = nan;
    res.AP_rise_time        = nan;
end

%% Find fall time of AP
Vm_postAP = Vm(res.AP_ind:res.fAHP_ind); 
amp_to_min = Vm_postAP(1) - Vm_postAP(end); 
t_postAP = dt * (0:length(Vm_postAP)-1);
t_postAP_interp = linspace(t_postAP(1), t_postAP(end), length(t_postAP)*fall_interp_factor);
Vm_postAP = interp1(t_postAP, Vm_postAP, t_postAP_interp);

shift_ind_post = res.AP_ind - 1;

fall_start_Vm           = fall_start_APfactor * amp_to_min + Vm_postAP(end);
fall_stop_Vm            = fall_stop_APfactor * amp_to_min + Vm_postAP(end);

fall_start_ind          = find_nearest_ind(Vm_postAP, fall_start_Vm);
fall_stop_ind           = find_nearest_ind(Vm_postAP, fall_stop_Vm);

res.AP_fall_start_time  = t_postAP_interp(fall_start_ind) + shift_ind_post*dt;
res.AP_fall_stop_time   = t_postAP_interp(fall_stop_ind) + shift_ind_post*dt;
res.AP_fall_start_Vm    = Vm_postAP(fall_start_ind);
res.AP_fall_stop_Vm     = Vm_postAP(fall_stop_ind);

res.AP_fall_time        = abs(res.AP_fall_stop_time - res.AP_fall_start_time);

%% Find half-width of AP
if ~isnan(res.AP_thres_ind)
    half_height_Vm = amp_to_min/2 + res.fAHP_Vm;
    res.AP_width_start_time = t_preAP_interp(find_nearest_ind(Vm_preAP, half_height_Vm))   + shift_ind_pre*dt;
    res.AP_width_stop_time  = t_postAP_interp(find_nearest_ind(Vm_postAP, half_height_Vm)) + shift_ind_post*dt;
    res.AP_width_Vm = half_height_Vm;
    res.AP_width = abs(res.AP_width_stop_time - res.AP_width_start_time);
else
    res.AP_width_start_time = nan;
    res.AP_width_stop_time  = nan;
    res.AP_width_Vm         = nan; 
    res.AP_width            = nan;
end
end
