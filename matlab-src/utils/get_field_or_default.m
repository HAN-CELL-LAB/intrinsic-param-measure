function val = get_field_or_default(s, field, def_val)
% Get field value if `field` exist, otherwise return default value `def_val`
% INPUT: 
%   s:          struct object 
%   field:      field name to look for 
%   def_val:    default value to return of field is not found
% OUTPUT:
%   val:        field value 

val = def_val;
if isfield(s, field)
    val = s.(field);
end

end