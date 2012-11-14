%% k = test_find_k('synthmod_01_01.csv', 3, 1 , 2, .5)


function k = test_find_k(data, i, j, r, spline_step)
data = csvread(data);              % read in csv file without 1st row, (labels)
data = data(2:end,:);   
time = data(:,1);                  % time vector is first column of data
samples = length(time);
nodes = data(:,2:end);             % redefine data without time column 
number_nodes = length(nodes(1,:)); % number of genes/species
splined_time = time(1):spline_step:time(end);
splined_nodes = zeros(length(splined_time),number_nodes);
for n = 1 : number_nodes
 splined_nodes(:,n) = interp1(time,nodes(:,n),splined_time,'spline'); % cubic spline interpolation
end
nodes = splined_nodes;
k = find_k(nodes(:,i), nodes(:,j), nodes(:,r), spline_step)
endfunction

 %k = find_k(nodes(:,1), nodes(:,2), nodes(:,3), 10)
