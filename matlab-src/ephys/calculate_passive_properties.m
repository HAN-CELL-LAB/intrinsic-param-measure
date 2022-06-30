function res = calculate_passive_properties(t, Vm, cell_info, stim_info, passive_params, dry_run)

if ~exist('dry_run', 'var')
    dry_run = false;
end

Vbase_props = calculate_Vbase(t, Vm, stim_info, passive_params.Vbase_calc_time_window, dry_run);
Rin_props = calculate_Rin(t, Vm, stim_info, passive_params.Rin_calc_time_window, dry_run);
Vrest_props = estimate_Vrest(Vbase_props.Vbase, cell_info.Ihold, Rin_props.Rin, dry_run);

res = mergefield_struct(Vbase_props, Rin_props, Vrest_props);

end

 
    