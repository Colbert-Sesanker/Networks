# A periodic distance meausre not mentioned in \ref{Cont per Osc} test
# Colbert Sesanker 2013

kernel_step = .01;  kernel_points = round(1/kernel_step); 
spline_step = .01;  spline_mesh = 0 : spline_step : kernel_points-spline_step; 
spline_points = kernel_points/spline_step;

window_size = 100; windows=spline_points/window_size;
window_vec = 1:windows;

for i=1:15
  # Cubic spline the frequency spectrum to better count dirac delta point masses
  periodic_spec = abs(fft(kernel(5,kernel_step,1)));
  splined_periodic_spec = interp1(0:kernel_points-1, periodic_spec, spline_mesh,'linear');  
  periodic_spec = splined_periodic_spec;

  random_spec =   abs(fft(kernel(2,kernel_step,1)));
  splined_random_spec = interp1(0:kernel_points-1, random_spec, spline_mesh,'linear');
  random_spec = splined_random_spec;

  smooth_spec =   abs(fft(kernel(3,kernel_step,1500)));
  splined_smooth_spec = interp1(0:kernel_points-1, smooth_spec, spline_mesh,'linear');
  splined_smooth_spec(splined_smooth_spec<0) = 0;
  smooth_spec = splined_smooth_spec;  
  
  nth_window = @(vec, n) vec(window_size*(n-1) + 1: window_size*n);
  integrated_windows = @(spec) arrayfun(
                                       @(n) spline_step*trapz(nth_window(spec,n)),
                                       window_vec
                                       );
  # Div approximates the sequence of dirac delta masses of frequency spectrum
  div = @(vec) arrayfun(@(x,y) max(x,y)/min(x,y), vec'(1:end-1), vec'(2:end));
  delta_masses = windows - 1;
  # Equation of the line passing through points: (1, min(deltas)) and (delta_masses, max(deltas))
  y = @(max, min) (@(x) x*(max-min)/(delta_masses-1) + (delta_masses*min - max)/(delta_masses-1));
  total_mass = @(max,min) quad(y(max,min), 1, windows-1);
 
  periodic_deltas =  sort(div(integrated_windows(periodic_spec))); 
  pd(:,i) =  periodic_deltas/periodic_deltas(end);
  max_p = periodic_deltas(end); min_p = periodic_deltas(1); 
  p_mass = total_mass(max_p, min_p); curved_mass_p = trapz(periodic_deltas); 
  periodic_measure(:,i) = p_mass/curved_mass_p;
 
  random_deltas =   sort(div(integrated_windows(random_spec)));
  rd(:,i) =  random_deltas/random_deltas(end);
  max_r = random_deltas(end); min_r = random_deltas(1); 
  r_mass = total_mass(max_r, min_r); curved_mass_r = trapz(random_deltas); 
  random_measure(i) = r_mass/curved_mass_r;
 
  smooth_deltas =   sort(div(integrated_windows(smooth_spec)));
  sd(:,i) =  smooth_deltas/smooth_deltas(end);
  max_s = smooth_deltas(end); min_s = smooth_deltas(1); 
  s_mass = total_mass(max_s, min_s); curved_mass_s = trapz(smooth_deltas); 
  smooth_measure(i) = s_mass/curved_mass_s;  
endfor

plot(periodic_measure - random_measure); plot.new();
plot(periodic_measure - smooth_measure)
