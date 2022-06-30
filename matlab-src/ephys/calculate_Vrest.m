function res = calculate_Vrest(t, Vm, stim_info, time_window)
if ~exist('time_window', 'var') 
    time_window = 10; 
end

time_end = get_field_or_default(stim_info, 'pulse_on_time', time_window); 

Vr = Vm(t>=time_end-time_window & t<=time_end);
res.Vrest = mean(Vr);
res.std_Vrest = std(Vr); 

end