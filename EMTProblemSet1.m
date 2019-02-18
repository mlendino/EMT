%% Michael Lendino EMT PSET 2
clc;
clear all;
%% Problem 1
%Consider a material with sigma=1.6, epsilon_r=4 and mu_r=1
%Find the frequencies where the loss tangent takes on values 10,2,1,1/2 and
%1/10
sigma = 1.6;
er = 4;
mu = 4*pi*(10^-7);
e0 = (1/(36*pi))*(10^-9);
e = er*e0;
%The loss tangent is defined as sigma/(omega*epsilon), where epsilon is
%epsilon in free space*er; LT = sigma/(omega*epsilon) =
%sigma/(omega*epsilon in free space*er) ==> omega = sigma/(LT*epsilon in free space*er)
%==> f = sigma/(LT*2pi*omega*epsilon in free space*er)
%uncomment value for loss tangent dependeing on what frequency you wish to
%calculate (so: 10,2,1,1/2,1/10)
%LT = 10 %==> f = .72GHz ==> alphadBexact = 557.24 ==>alphaapprox = 585.76
%LT = 2 %==> f = 3.6GHz ==> alphadBexact = 1.03e3 ==>alphaapprox = 1.31e3
%LT = 1 %==> f = 7.2GHz ==> alphadBexact = 1.19e3 ==>alphaapprox N/A 
%LT = 1/2 %==> f = 1.44e10 Hz ==> alphadBexact = 1.27e3 ==>alphaapprox = 1.31e3
LT = 1/10 %==> f = 7.2e10Hz ==> alphadBexact = 1.31e3 ==>alphaapprox = 1.31e3
f = sigma/(2*pi*LT*e);
w = 2*pi*f;
%at each frequency also compute the exact attenuation constant alpha in
%dB/m and the exact wavelength; the exact formula is gamma =
%j*w*sqrt(mu*e)*sqrt[1+ (sigma/j*w*e)] = alpha + jBeta
gamma = j*w*sqrt(mu*e)*sqrt(1+(sigma/(j*w*e)));
alpha = real(gamma);
alphadB = (20*log10(exp(1)))*alpha;
%also for the cases where the loss tangent>1 compute the results using the
%good conductor approximation: alpha = sqrt(pi*f*mu*sigma)
alphagoodc = sqrt(pi*f*mu*sigma);
%alpha = beta for a good conductor
alphagoodcdB = (20*log10(exp(1)))*alphagoodc;
%and for loss tangent<1, compute the results using the good dielectric
%approximation: alpha = 0.5*sigma*sqrt(mu/epsilon) and beta = w*sqrt(mu*e)*(1+((1/8)*((sigma/(w*e))^2)))
alphagoodd = 0.5*sigma*sqrt(mu/e);
alphagoodddB = (20*log10(exp(1)))*alphagoodd;
betagoodd = w*sqrt(mu*e)*(1+((1/8)*((sigma/(w*e))^2)))
%present all results in a table
%we get the exact wavelength from doing lambda = 2*pi/(beta)
beta = imag(gamma);
lambda = (2*pi)/beta;
lambaapproxgoodc = (2*pi)/alphagoodc;
lambaapproxgooddielectric = (2*pi)/betagoodd

LossTangent = [10;2;1;1/2;1/10];
FrequencyInGHz = [.72;3.6;7.2;14.4;72];
WavelengthInMetersEXACT = [0.0886;.0328;.0190;.0101;.0021];
%dBperm means dB/m
AlphaExactdBperm = [557.24;1.03e3;1.19e3;1.27e3;1.31e3];
AlphaApproxdBperm = [587.76;1.31e3;1.19e3 ;1.31e3;1.31e3];
WavelengthApproxInMeters = [.0932;.0417;.0190;.0101;.0021];

T = table(LossTangent,FrequencyInGHz,WavelengthInMetersEXACT,WavelengthApproxInMeters,AlphaExactdBperm,AlphaApproxdBperm)











