function moment = spectral_moment(f,S,f_1,f_2,p)
% Compute the p^th spectral density moment
% withing the frequency band [f_1 f_2]
% Inputs: 
%    f     frequency vector in Hz associated with S 
%    S     variance density spectrum ([unit x]^2/Hz, where x is the initial timeseries)
%    f_1   lower limit frequency (Hz)
%    f_2   upper limit frequency (Hz)
%    p     order of the moment 
%          (compute the zeroth-order moment is p is not specified)
% Output:  
%    moment p^th order spectral moment ([unit x]^2 Hz^p)

% if p not specified, p=0 by default
if nargin < 5,  
   p = 0;
end;

% Determination of the elements of f contained between f_1 and f_2
index = (f<=f_2)&(f>=f_1);

% Computation of the p^th order moment
m = (f(index).^p).*S(index);
moment =  trapz(f(index),m);


   