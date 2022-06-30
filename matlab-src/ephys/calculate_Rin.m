function res = calculate_Rin(t, Vm, stim_info, time_window, dry_run)
% V [mV]
% I [pA]
% R [MOhm]

res = struct;
res.units.V = 'mV';
res.units.I = 'pA';
res.units.R = 'MOhm';
res.dV_forRin = nan;
res.std_VforRin_before = nan;
res.std_VforRin_after = nan;
res.Rin = nan;

if ~exist('dry_run', 'var')
    dry_run = false;
end

if dry_run
    return;
end


if ~exist('time_window', 'var') 
    time_window = 10; 
end
necessary_fields = {'hyperpol_on_time', 'hyperpol_off_time'};
if ~all(cellfun(@(x) isfield(stim_info, x), necessary_fields))
    error('Need on and off time of hyperpolarizing currents'); 
end

hyperpol_on_time  = stim_info.hyperpol_on_time;
hyperpol_off_time = stim_info.hyperpol_off_time;

hyperpol_Iinj_amp = get_field_or_default(stim_info, 'hyperpol_Iinj_amp', nan); 

Vm_before = Vm(t>=hyperpol_on_time-time_window & t<=hyperpol_on_time);
Vm_after  = Vm(t>=hyperpol_off_time-time_window & t<=hyperpol_off_time);


res.dV = abs(mean(Vm_after) - mean(Vm_before));
res.std_V_before = std(Vm_before); 
res.std_V_after = std(Vm_after); 

res.Rin = (res.dV * 1e-3) / (hyperpol_Iinj_amp * 1e-12) /1e6;

end