function v = getfield_structarray(s, field, uni_out)
if ~exist('uni_out', 'var'), uni_out = 1; end
try 
    v = arrayfun(@(x) x.(field), s, 'uni', uni_out);
catch
    warning('The output is not uniform, hence will turn ''Uni'' to 0');    
    v = arrayfun(@(x) x.(field), s, 'uni', 0);
end
end