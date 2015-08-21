function [ freq, spectral_power ] = Fourier_Transform( time_domain_signal, time_step )
%%%
%
% Performs an FFT on the time domain signal.
%
% Ryan Skinner, August 2015
%
%%%

% Clean up the signal. A single NaN makes fft return all NaNs.
time_domain_signal = time_domain_signal(~isnan(time_domain_signal));

Fs = 1/time_step;    % Sampling frequency
Fnyquist = Fs/2;    % Nyquist frequency

L = length(time_domain_signal);   % Length of signal (steps of period T=1/Fs)
nFFT = 2^nextpow2(L);           % For efficiency and to use all data

% Compute spectral power, which is the imaginary number's magnitude
% Also only use the first half of the transform; the other is mirrored
spectral_power = abs(fft(time_domain_signal, nFFT)).^2;
spectral_power = spectral_power(1:nFFT/2);

% Frequency is evenly spaced from 0 to the Nyquist frequency.
freq = Fnyquist * linspace(0, 1, nFFT/2);

end