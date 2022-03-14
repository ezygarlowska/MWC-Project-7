function [S f edf conf95Interval] = VarianceDensitySpectrum(x,nfft,Fs)
% compute variance spectral density of series x 
%
% computational steps:
%   1. entire series is detrended (second-order polynomial)
%   2. series is divided into blocks with length nfft, using 50% overlap
%   3. each block is tapered with Hamming window of length nfft
%   4. fft, conversion from one- to two-sided estimate
%   5. 95% confidence interval (provided chi2inv.m is available)
%
% INPUT
%   x, timeseries [n x 1]
%   nfft, block length (even number)
%   Fs, sampling frequency (Hz)
%
% OUTPUT
%   S, variance spectral density, (unit x)^2/Hz
%   f, frequency axis (Hz)
%   edf, effective degree of freedom
%   conf95Interval, [1 x 2], lower-upper 95% multiplication factor for S
%
% v1, Gerben Ruessink, 30-10-2010.
% Based on Bendat and Piersol (2000)

% some initial preparations
x = x(:);
n = size(x,1);
nfft = nfft - rem(nfft,2);

% detrend the series
tDummy = (1:n)';
p = polyfit(tDummy,x,2); 
x = x - polyval(p,tDummy);

% prepare the blocks
overlap = 0.5;
nOverlap = overlap*nfft;
nBlocks = fix((n-nOverlap)/(nfft-nOverlap));
xb = NaN(nfft,nBlocks);
id = 1:nfft;
for b = 1:nBlocks
    xb(:,b) = x(id);
    id = id + (nfft-nOverlap);
end;
clear id

% taper each block
ww = window(@hamming,nfft);
normFactor = mean(ww.^2);
for b = 1:nBlocks
    xb(:,b) = xb(:,b).*ww / sqrt(normFactor); % to ensure that var(xb) before windowing == var(xb) after windowing
end;

% frequency axis
f = 0:nfft/2;
fId = f + 1;
df = Fs/nfft;
f = f'*df;

% fft
S = fft(xb,nfft,1)/nfft;
S = S(fId,:);                % only one side needed
S = mean(abs(S).^2,2)/df;    % average over blocks and normalize
S(2:end-1) = 2*S(2:end-1);   % from one- to two-sided

% confidence statistics
edf = round(nBlocks*2.52164); % Based on Priestley (1981), Table 6.2, page 467
if exist('chi2inv','file'),
    alpha = 0.05;
    conf95Interval(1) = edf/chi2inv(1-alpha/2,edf);
    conf95Interval(2) = edf/chi2inv(alpha/2,edf);
else
    conf95Interval(1:2) = NaN;
end;