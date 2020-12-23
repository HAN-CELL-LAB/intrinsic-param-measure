function [x_nearest, ind] = find_nearest(v, x, return_first)
if ~exist('return_first', 'var'), return_first = 'val'; end
[~, ind] = min(abs(v-x));
x_nearest = v(ind);
if strcmpi(return_first, 'val') % do nothing 
elseif strcmpi(return_first, 'ind')
    [x_nearest, ind] = swap(x_nearest, ind); 
else 
    error('The parameter "return_first" can only be "val" or "ind" to indicate what to return first');
end
end