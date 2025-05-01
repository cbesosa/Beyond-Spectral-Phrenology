function [frequency, power, all_in_one] = psd1(real_signal, signal_frequency)


nt = length(real_signal);
[power, frequency] = cpsd(real_signal, real_signal,hamming(nt),0,nt,signal_frequency);

% power is withrespect to the windowed signal.

% power = power./2.517;

% check consistency of power:
%{ 
X = real_signal.*hamming(nt);   
[power(1),  mean(X).^2./(length(X)/signal_frequency),   power(1)/(mean(X).^2)]
[sum(power(2:end)),  var(X),   sum(power(2:end))/var(X)]  
%}

all_in_one.f = frequency;
all_in_one.p = power;
all_in_one.w = length(real_signal)/signal_frequency;  %  in unit of time.


end