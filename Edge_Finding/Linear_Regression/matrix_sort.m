%% Colbert Sesanker 10/24/2012
%% Sorts all elements in matrix into vector and returns
%% Indicies of sorted elements in original matrix

function [A_sorted, index] = matrix_sort(A, num_elts)
ll = [];
cols = size(A,2);
for i = 1 : size(A,1)
  ll = [ll, A(i,:)];
end
[sorted, idx] = sort(ll);
A_sorted = sorted(1:num_elts);
index = [];
for i = 1:length(idx)
  row = ceil(idx(i)/cols);
  column = idx(i) - (row-1)*cols;   
  index = [index; [row,column]];
end
index = index(1:num_elts,:);
endfunction


