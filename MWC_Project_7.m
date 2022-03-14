clear; clc; close all;

%% Chapter 3 - Spectral Analysis

%% 3.1 Computation of the density spectrum 

%% Load the data 

%Loading the files LowTide.txt, midTide.txt and highTide.txt
ucData = load('ucData.mat');
condition1=ucData.condition1;
condition2=ucData.condition2;
condition3=ucData.condition3;
%detrended time-series of free- surface elevation (m)
eta1 = condition1.data(:,1);
eta2 = condition2.data(:,1);
eta3 = condition3.data(:,1);
%cross-shore near-bed velocity (m/s) (Note that positive values of u
%indicate onshore-directed velocities.)
u1 = condition1.data(:,2);
u2 = condition2.data(:,2);
u3 = condition3.data(:,2);
%near-bed sediment concentration (kg/m3)
conc1 = condition1.data(:,3);
conc2 = condition2.data(:,3);
conc3 = condition3.data(:,3);
%water depth
h1 = condition1.h;
h2 = condition2.h;
h3 = condition3.h;

%% Compute the variance density spectrum 

% Computing the variance density spectra of the free-surface elevation for the three conditions.
% Number of blocks is equal 15, so that we have 30 degrees of freedom
x= [eta1 eta2 eta3];
Fs=condition1.Fs; %it is the same for each of the conditions
figure;
for i=1:3
    subplot(3,1,i);
    [S f edf conf95Interval] = VarianceDensitySpectrum(x(:,i),0.5^3*length(x(:,i)),Fs);
    plot(f,S);
    grid on;
    %ylim([0 9]);
    if i==1 
        title('Variance density spectra of the free-surface elevation');
    end 
    %ylabel('E [m^2/Hz]','FontWeight','bold');
    legend('Spectrum (N=15)');
end
%xlabel('Frequency [Hz]','FontWeight','bold');

%add labels and confidence intervals, remember that we did something wrong
%in second practical, so go to the feedback and fix conf95 thingy

%% spectral wave height

%For each condition, compute the spectral wave height Hm0 for the sea-swell (0.05 Hz < f < fN, with fN the Nyquist frequency) and infragravity wave (0.005 Hz < f < 0.05 Hz) frequencies.

%sea-swell 0.05 Hz < f < fN
fmin=0.05;
fmax=1/(2*0.5);
    for i=1:3
        [S f edf conf95Interval] = VarianceDensitySpectrum(x(:,i),0.5^3*length(x(:,i)),Fs);
        
        m0(i) = spectral_moment(f,S,fmin,fmax);
        H_m01 = 4*sqrt(m0);   %So now H_m0 = 4*sqrt(m_0)

    end

%infragravity waves 0.005 Hz < f < 0.05 Hz
fmin=0.005;
fmax=0.05;
    for i=1:3
        [S f edf conf95Interval] = VarianceDensitySpectrum(x(:,i),0.5^3*length(x(:,i)),Fs);
        
        m0(i) = spectral_moment(f,S,fmin,fmax);
        H_m02 = 4*sqrt(m0); 

    end

H_m0= [H_m01 ;H_m02];
display(H_m0);

%% relative wave height

%Compute the relative wave height, defined as Hm0,ss/h, and the ratio Hm0,inf /Hm0,ss for each case.

h=[h1 h2 h3];

for i=1:3
    rwh(i)=H_m0(1,i)/h(i); %relative wave height Hm0,ss/h
    ratio(i)=H_m0(2,i)/H_m0(1,i); %ratio Hm0,inf /Hm0,ss
end 
display(rwh);
display(ratio);

%we should maybe make a nice plot with all the characteristics or smth

%% skewness and asymmetry

%Compute the skewness and asymmetry of the sea-swell waves for each case.

flow=0.05;
fhigh=1/(2*0.5);
for i=1:3
    dataFilt (:,i)= fft_filter(x(:,i), Fs,flow , fhigh);
end 


