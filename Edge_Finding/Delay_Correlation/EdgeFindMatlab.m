%%% First read csv file into matrix and extract nodes. 
%%% Chop off labels & time, Edge Finding Algorithm based on Delays
%%% Colbert Sesanker 09/2012

data = csvread('synthmod_01_01.csv', 1); % read csv excluding first row, the labels 
time = data(:,1);                 % time vector is first column of data
time_pts = length(time);          % number of time points
nodes = data(:,2:end);             % redefine data without time column                                 
number_nodes = length(nodes(1,:)); % take length of arbitrary column



num_pairs = nchoosek(number_nodes - 1, 2); % pairs between every node other than the ith examined node
num_functions = 4;                         % update as functions are added
time_delays = 3;
ed = 40;  % Expected Delay

%%% 'correlations' is matrix of matricies. The ith matrix offers explanations of node i
%%% in terms of correlations to functions (max, min ...) evaluated pairwise on
%%% the other nodes

correlations = zeros(num_functions, num_pairs, number_nodes, time_delays); 
single_node_correl = zeros(number_nodes, number_nodes - 1, time_delays); % compute correlations with the other individual nodes
for i = 1 : number_nodes
  for t = ed : ed : 3*ed
    node = nodes(1+t/10:end,i); % the node being examined for explanation (boundary conditions are actually periodic)
    l = 1;             % l is a counter that counts over the number of pairs
    for j = 1 : number_nodes    
      if (j ~= i)         
          % correlate node i with jth node (not equal i). equivalent to covariance-matrix - eye
          single_node_correl(i,j,t/ed) = cor(node, nodes(1:end-t/10,j));
        for k = 1 : number_nodes
          if (k ~= i && k > j)        % nodes cannot explain themselves and pairs have no order 
            p = corr(max(nodes(1:end-t/10,j), nodes(1:end-t/10,k)), node);
            correlations(1,l,i,t/ed) = cor(max(nodes(1:end-t/10,j), nodes(1:end-t/10,k)), node); % correlate node i with all pairs != i
            correlations(2,l,i,t/ed) = cor(min(nodes(1:end-t/10,j), nodes(1:end-t/10,k)), node); % using each pairwise function
            correlations(3,l,i,t/ed) = cor(nodes(1:end-t/10,j) + nodes(1:end-t/10,k), node);
            correlations(4,l,i,t/ed) = cor(nodes(1:end-t/10,j) .* nodes(1:end-t/10,k), node);  
            l = l + 1;
          else
           continue
          end
        end
      else
        continue
      end
    end
  end
end
