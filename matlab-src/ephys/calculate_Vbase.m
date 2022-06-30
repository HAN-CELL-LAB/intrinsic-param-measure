function res = calculate_Vbase(t, Vm, stim_info, time_window, dry_run)

res = struct;
res.Vbase = nan;
res.std_Vbase = nan;

if ~exist('dry_run', 'var')
    dry_run = false;
end

if dry_run
    return;
end

if ~exist('time_window', 'var') 
    time_window = 10; 
end

time_end = get_field_or_default(stim_info, 'pulse_on_time', time_window); 

Vb = Vm(t>=time_end-time_window & t<=time_end);
res.Vbase = mean(Vb);
res.std_Vbase = std(Vb); 

end