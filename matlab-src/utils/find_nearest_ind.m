function ind = find_nearest_ind(v, x)
[~, ind] = min(abs(v-x));
end