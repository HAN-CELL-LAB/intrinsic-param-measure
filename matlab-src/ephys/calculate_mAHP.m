function res = calculate_mAHP(Vm, ind_start, dt, params)
% TODO: edge case when spike happens at the beginning 
% Vm [mV]
% ind_AP [no unit]
% dt [ms]
% params: struct
%     max_tpost [ms] [def=300] : max time to look forward
%     find_opt [string] [def=min]
%         - 'raw': get min from raw 
%         - 'smooth': smoothed then min
%         - 'doubexp': fitted with double exp then min
%         - 'lowpass': low pass then min 
%         - 'all': return all 
%     save_processed [bool] [def=false]: whether to save the processed
%     (smoothed, fitted or filtered) data as well


allowed_opts = {'raw', 'smooth', 'doubexp', 'lowpass', 'all'}; 
find_opt = params.find_opt; 

max_tpost = get_field_or_default(params, 'max_tpost', 300); 
save_processed = get_field_or_default(params, 'save_processed', false); 

incr_post = ceil(max_tpost/dt); 
Vm_post = Vm(ind_start:ind_start+incr_post); 

method_info = struct; 

switch lower(find_opt)
    case 'raw'
        method_info.description = 'raw -> min'; 
    case 'smooth'
        smooth_method = get_field_or_default(params, 'smooth_method', 'movmean');
        smooth_tau    = get_field_or_default(params, 'smooth_tau', 5); % ms
        
        Vm_post = smoothdata(Vm_post, smooth_method,  ceil(smooth_tau/dt));
        method_info.description = sprintf('smooth (%s,%.1f ms) -> min', smooth_method, smooth_tau); 
    case 'doubexp' 
        doubexp_lowbound  = get_field_or_default(params, 'doubexp_lowbound', [-50,-100,0,0,-90]);
        doubexp_uppbound  = get_field_or_default(params, 'doubexp_uppbound', [100,50,150,500,-40]);
        doubexp_startpnt  = get_field_or_default(params, 'doubexp_startpnt', [10,-10,50,100,-60]);
        
        fit_opt_obj = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',doubexp_lowbound,'Upper',doubexp_uppbound,...
            'StartPoint',doubexp_startpnt);
        doubexp_eq = 'a1*exp(-x/tau1) + a2*exp(-x/tau2) + y_c'; 
        fit_type_obj = fittype(doubexp_eq, 'options', fit_opt_obj);
        
        t = (0:length(Vm_post)-1)'*dt; 
        
        fit_fun = fit(t,Vm_post,fit_type_obj);
        method_info.fit_params = struct('a1', fit_fun.a1, 'a2', fit_fun.a2, ...
            'tau1', fit_fun.tau1, 'tau2', fit_fun.tau2, 'y_c', fit_fun.y_c); 
        Vm_post = fit_fun(t);
        
        method_info.description = 'double exp fit -> min'; 
    case 'lowpass'
        lowpass_freq = get_field_or_default(params, 'lowpass_freq', 50); 
        Fs = 1e3/dt;
        
        [b, a] = butter(2, lowpass_freq/(Fs/2), 'low');
        Vm_post = filtfilt(b,a,Vm_post);
        method_info.description = sprintf('lowpass %.1f Hz -> min', lowpass_freq);
    case 'all'
        for i = 1:(length(allowed_opts)-1)
            sub_opt = allowed_opts{i}; 
            params.find_opt = sub_opt; 
            res.(sub_opt) = calculate_mAHP(Vm, ind_start, dt, params); 
        end
        return; 
    otherwise
        error('"%s" is not an allowed option. Only one of these are allowed \n\t%s\n', ...
            find_opt, sprintf('"%s"  ', allowed_opts{:}));
end

[res.mAHP_amp, res.mAHP_ind] = min(Vm_post);
res.mAHP_ind = res.mAHP_ind + ind_start - 1;
res.mAHP_time = (res.mAHP_ind-1) * dt; 
res.method_info = method_info; 

res.processed_Vm = nan; 
res.processed_t = nan; 
if save_processed
    res.processed_Vm = Vm_post; 
    res.processed_t = ((ind_start-1) + (0:length(Vm_post)-1))*dt;
end

end