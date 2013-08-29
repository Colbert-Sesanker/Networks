# Colbert Sesanker 3/2013
# Tests periodic distance (pd) measure for robustness
# Should be invariant to frequency (freq), amplitude (amp),  
# phase, translation (dc_comp) and gaussian noise
# Should vary continously and monotonically with freq
# if desired can be fed a gaussian proceess

function score, freq_sens_curve = periodic_distance_test()

freq_sens = @(frequency) pd_sensitivity(freq=frequency, amp=1, phase=0, 
                                        dc_comp=0, noise_scale=0)
freq_sens_curve = arrayfun(freq_sens, [1: .1: 10]);
figure(1);
plot([0: .1: 9], freq_sens_curve);
title("Frequency Sensitivity");
xlabel("Frequency Difference between Signals (Hz)");
ylabel("Periodic Distance between Signals");

amp_sens = @(amplitude) pd_sensitivity(freq=1, amp=amplitude, phase=0, 
                                       dc_comp=0, noise_scale=0);
amp_sens_curve = arrayfun(amp_sens, [1: 10: 1000]);
figure(2);
plot([0: 10: 999], amp_sens_curve);
title("Amplitude Sensitivity");
xlabel("Amplitude Difference between Signals");
ylabel("Periodic Distance between Signals");

phase_sens = @(phase_) pd_sensitivity(freq=1, amp=1, phase=phase_, 
                                      dc_comp=0, noise_scale=0);;
phase_sens_curve = arrayfun(phase_sens, [0: .1: pi]);
figure(3);
plot([0: .1: pi], phase_sens_curve);
title("Phase Sensitivity");
xlabel("Phase Difference between Signals");
ylabel("Periodic Distance between Signals");

dc_sens = @(dc) pd_sensitivity(freq=1, amp=1, phase=0, 
                               dc_comp=dc, noise_scale=0);;
dc_sens_curve = arrayfun(dc_sens, [0: 10: 1000]);
figure(4);
plot([0: 10: 1000], dc_sens_curve);
title("DC offset Sensitivity");
xlabel("DC offset difference between Signals");
ylabel("Periodic Distance between Signals");

noise_sens = @(scale) pd_sensitivity(freq=1, amp=1, phase=0, 
                                     dc_comp=0, noise_scale=scale);;
noise_sens_curve = arrayfun(noise_sens, [0: .1: 10]);
figure(5);
plot([0: .1: 10], noise_sens_curve);
title("Gaussian Noise Sensitivity");
xlabel("Noise difference between Signals");
ylabel("Periodic Distance between Signals");

# Compute sum of changes to all invariants. Should be close to zero
score = sum([amp_sens_curve phase_sens_curve dc_sens_curve noise_sens_curve])
endfunction

function sens = pd_sensitivity(freq=1, amp=1, phase=0, 
                               dc_comp=0, noise_scale=0)
time = [-10*pi: .01: 10*pi];
noise = randn(size(time));

sens = periodic_distance(
                       amp*cos(2*pi*freq*time + phase) + dc_comp + noise_scale*noise,
                       cos(2*pi*time)
                       );
endfunction

function sens = pd_sens_complex(freq=1, amp=1, phase=0, 
                                dc_comp=0, noise_scale=0)
time = [-10*pi: .01: 10*pi];
noise = randn(size(time));

sens = periodic_distance(

amp*    cos(2*pi*    freq*time + phase       ) + dc_comp  + ...
amp* .5*cos(2*pi*.5* freq*time + phase + pi/2) + dc_comp  + ...    
amp*.25*cos(2*pi*.25*freq*time + phase + pi/4) + dc_comp  ,                                             
                                                                               
cos(2*pi*time) + .5*cos(2*pi*.5*time + pi/2) + .25*cos(2*pi*.25*time + pi/4)
                       );
endfunction

