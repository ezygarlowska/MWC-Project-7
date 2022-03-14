function dataFilt = fft_filter(data, fs, flow, fhigh)

% fft_filter is a spectral filter which basically sets the undesired
% frequencies to 0.
% flow and fhigh determine if the filter is a low-pass, high-pass or band-pass filter. 
% if flow = 0 :     low pass
% if fhigh = fNyq : high pass
% else band-pass
%
% INPUT
%   data      initial timeseries 
%   fs        sampling frequency of the time-series (Hz)
%   flow      lowest frequency desired (Hz)
%   fhigh     highest frequency desired (Hz)
% OUTPUT
%   dataFilt  filtered timeseries 



data = data(:);
[Nt, Nx] = size(data);
df = fs/Nt;

% fft (fast fourier transform) of data
Y = fft(data);

% define the frequency axis
odd = (round(Nt/2)-(Nt/2) == 0.5);
if odd
   f = [0:(Nt-1)/2];
   f = [f , -(Nt-1)/2+f(1:(Nt-1)/2)];
else
   f = [0: Nt/2];
   f = [f,-Nt/2+f(2:(Nt/2))];
end;
f = f(:)*df;

% look for the frequency which should be removed 
t = ((abs(f)< flow) | (abs(f) > fhigh));

% set them to zeros
Y(t) = 0;

% back to time domain (inverse fourier transform)
dataFilt = real(ifft(Y));

