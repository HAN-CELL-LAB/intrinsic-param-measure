function ind = find_start_of_sustained_on(v, min_pts) 
% Find the index of sustained on for `min_pts` points of boolean vector `v`
% EXAMPLE: 
%   with `min_pts=3` and
%   >>> v = [0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1]
%   Then `ind_start` would be 9
% INPUT:
%   v:          boolean vector 
%   min_pts:    minimum # of points required to be `true`
% OUTPUT:
%   ind:        the first index of sustained on (`true`) for `min_pts` 
% NOTE: inspired by https://stackoverflow.com/a/29329330

x = filter(ones(min_pts,1),1,v); 
ind = find(x == min_pts, 1) - (min_pts-1);

end