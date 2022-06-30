function res = estimate_Vrest(Vbase, Ihold, Rin, dry_run)
% Vbase [mV]
% Ihold [pA]
% Rin [MOhm]

res = struct;
res.Vrest_1 = nan;
res.Vrest_2 = nan;

if ~exist('dry_run', 'var')
    dry_run = false;
end

if dry_run
    return;
end


Vbase = Vbase * 1e-3;
Ihold = Ihold * 1e-12;
Rin = Rin * 1e6; 

Vr_1 = Vbase - Ihold * Rin;

Ix = Vbase/Rin;
Vr_2 = Vbase*Ix/(Ix + Ihold);  

res.Vrest_1 = Vr_1 * 1e3; 
res.Vrest_2 = Vr_2 * 1e3; 

end