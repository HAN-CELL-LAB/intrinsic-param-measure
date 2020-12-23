function v = return_nan_if_empty(v)
if isempty(v) 
    v = nan; 
end
end