%% Colbert Sesanker, Jan 2013
%% This example plots samples from a Gaussian Process for 5 kernels:
%% load function into workspace. 
%% kernel(3, .005, 100) will plot a sample with guassian kernel, grid step = .005 and c = 100 

function GP_sample = kernel(type, grid_step, c, plt=false)  
kernel =  {@(x,y)c*x'*y, % linear type = 1
           @(x,y)c*min(x,y), % brownian motion, type = 2 
           @(x,y)exp(-c*(x-y)'*(x-y)),   % gaussian, type = 3
           @(x,y)exp(-c*abs(x-y)),       % laplace,  type = 4
           @(x,y)exp(-c*sin(20*pi*x-y))}; % periodic, type = 5

% Create Kernel
x = 0:grid_step:1-grid_step;
n = length(x);
K = zeros(n,n);
for i = 1:n
    for j=1:n
        K(i,j) = kernel{type}(x(i), x(j));
    end
end

% Plot
u = randn(n,1);     % u ~ N(0, I)
[A, S, _] = svd(K); % ~ = A' and is not used
GP_sample = A*sqrt(S)*u;    % z ~ N(0, K) (by the affine property)
if plt == true,
  figure(1); 
  plot(x,z,'.-');
  axis([0,1,-2,2]);
end
endfunction


