function res = structarray_to_struct(structobj, uni_out)
if ~exist('uni_out', 'var'), uni_out = 1; end 
fields = fieldnames(structobj); 
res = struct(); 
for i = 1:length(fields)
    field_name = fields{i}; 
    res.(field_name) = getfield_structarray(structobj, field_name, uni_out);
end
end