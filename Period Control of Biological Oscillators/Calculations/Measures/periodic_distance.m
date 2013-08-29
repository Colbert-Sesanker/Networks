# This functiona computes the periodic distance between signals and is interfaced to python via 
# oct2py. Colbert Sesanker 3/2013

# Compares the k largest peaks in the period
# Signal_1 is base truth and fixed
function distance = periodic_distance(signal_1, signal_2, k=3, time_interval=[0,1000], amp=true)

N = length(signal_1);
assert(length(signal_1) == length(signal_2));

# Coerce input signals into row vectors
signal_1 = make_row_vector(signal_1); signal_2 = make_row_vector(signal_2);

# Normalize the signal sum of squares to 1 
normalize = @(signal) (signal - mean(signal)) ./ ...
             (std(signal) * sqrt(N - 1)) ;

signal_1 = normalize(signal_1); signal_2 = normalize(signal_2);

# Remove NaN and inf artifacts from normalization procedure
signal_1 = clean(signal_1);  signal_2 = clean(signal_2);

# Catch if the length of frequency_indexes is less than k
try
    # Find Dominant Frequencies and indicies of Signals 1, 2
    [dom_freqs_1, idx_exp] = dominant_frequencies(signal_1, k)
    [dom_freqs_2, idx_calc] = dominant_frequencies(signal_2, k)    
catch
    printf ("---WARNING---: Spectrum contains less than k peaks:...
            %s\n", k, lasterr)  
end
# Frequency distance in index units
delta_omega = 500;
vec_max  = @(vec1, vec2) arrayfun(@max, vec1, vec2);
# Control Switches
amp_diff = @(sig_1, sig_2) abs((max(sig_1) - min(sig_1)) - (max(sig_2) - min(sig_2))) / ... 
                               (max(sig_1) - min(sig_1))

# Allowed to vary by level percent of the native amplitude                              
amp_diff_switch = @(percent, n)  (amp_diff(signal_1(round(N/2):end), signal_2(round(N/2):end)) + (1 - percent))^n

freq_distance = @(a,b,c,d) sum( (abs(a-b) + .35).^inf + 
                                 exp(vec_max(a, b) .* abs(c - d - delta_omega))
                              );

distance = freq_distance(dom_freqs_1, dom_freqs_2, idx_calc, idx_exp) + amp_diff_switch(.5, inf)
endfunction

function signal_row = make_row_vector(signal)
if size(signal)(1) ~= 1
    signal_row = signal' ;
else
    signal_row = signal;
endif
endfunction

function cleaned = clean(signal)
signal(isnan(signal)) = 0; signal(isinf(signal)) = 0;
cleaned = signal;
endfunction

function power_spec = power_spectrum(signal)
num_points   = length(signal);
squared_spec = abs(fft(signal)).^ 2 / num_points;
power_spec   = 2*squared_spec(1: ceil(num_points / 2));
endfunction

function [dominant_freqs, dominant_freq_indexes] ...
          = dominant_frequencies(signal, k)

power_spec = power_spectrum(signal);

# Sort frequencies of power spectrum in decending order
# then take the k largest frequencies
[frequencies, ... 
 frequency_indexes] = sort(power_spec, 'descend');

dominant_freq_indexes = frequency_indexes(1:k);
dominant_freqs = frequencies(1:k);
endfunction









