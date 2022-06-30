function s = mergefield_struct(varargin)

if ~all(cellfun(@(x) isstruct(x) && length(x) == 1, varargin))
    error('All of inputs need to be a struct, and not arrays!');
end

all_fields = cellfun(@fieldnames, varargin, 'uni', 0); 
all_fields = vertcat(all_fields{:}); 
if length(all_fields) ~= length(unique(all_fields))
    error('There are overlapping fields across some structs!');
end

all_values = cellfun(@struct2cell, varargin, 'uni', 0); 
all_values = vertcat(all_values{:}); 

s = struct; 
for i = 1:length(all_fields) 
    s.(all_fields{i}) = all_values{i};
end
end

