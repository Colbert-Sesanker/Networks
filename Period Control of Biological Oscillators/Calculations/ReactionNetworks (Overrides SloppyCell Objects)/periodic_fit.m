function fit = periodic_fit(signal)
num_points = length(signal); 
if  num_points < 10000,
   signal_mesh = 1:num_points;
   spline_points = 10000;
   spline_step = (num_points-1) / (spline_points-1); 
   spline_mesh = 1: spline_step : num_points;
   # Cubic spline the time signal to better count dirac delta point masses
   spec = abs(fft(signal));
   # Cubic spline the time signal to better count dirac delta point masses
   splined_spec = interp1(signal_mesh, spec, spline_mesh, 'linear');
   spec = splined_spec;
   num_points = spline_points;
else
    nearest_hundred = num_points - mod(num_points, 100);
    signal = signal(1:nearest_hundred);
    spec = abs(fft(signal));
    num_points = nearest_hundred;
endif
# returns   
window_size = 100; windows = num_points / window_size;
window_vec = 1:windows;
nth_window = @(vec, n) vec(window_size*(n-1) + 1: window_size*n);

# Integrates freq spectrum window by window returning a vector of integrals
# The window size is small in an attempt to catch the 
# dirac delta masses of the spectrum
integrated_windows = @(spec) arrayfun(
                                     @(n) spline_step*trapz(nth_window(spec,n)),
                                     window_vec
                                     );
# Div amplifies the boundaries where the dirac delta masses spike
div = @(vec) arrayfun(@(x,y) max(x,y)/min(x,y), vec'(1:end-1), vec'(2:end));

spec_deltas = div(integrated_windows(spec));
sorted_deltas = sort(spec_deltas);
max_delta = sorted_deltas(end);
spec_deltas =  sorted_deltas/max_delta;
fit = spec_deltas;
endfunction







