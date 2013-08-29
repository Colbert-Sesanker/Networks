# Colbert Sesanker 3/2013
# calculates the difference in oscillation between signal_1, signal_2
# This is called in python using the  oct2py interface

function distance = oscillate_distance(signal_1, signal_2, k=5, time_interval=[0,1000])

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
    [dom_freqs_1, idx_exp] = dominant_frequencies(signal_1, k);
    [dom_freqs_2, idx_calc] = dominant_frequencies(signal_2, k);  
catch
    printf ("---WARNING---: Spectrum contains less than k peaks:...
            %s\n", k, lasterr)  
end
# Scale depends on 'level', relative to initial cost of 1: 
# exp(scale*level) is the cost jump from zero
# Probability of acceptance is exp(-distance)
level = .3;
scale = 10;
hill     =     @(delay, n, x) x.^n ./ (delay + x.^n);
sum_diff =     @(a, b)        abs(a - b);
switch_diff =  @(a, b, level) (sum_diff(a,b) + (1-level))
osc_distance = @(a,b,c,d)     sum(exp(scale*sum_diff(a,b) .* hill(1, 100, switch_diff(a, b, level))));
distance = osc_distance(dom_freqs_1, dom_freqs_2, idx_calc, idx_exp);
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









