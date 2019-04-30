%% Michael Lendino EMT PSET 7
clc;
clear all;
close all;
%% Question 1
%graph f(f/fc) vs f/fc 
freq = linspace(1,8,1000);
%from part a
dffc = freq.^-3./sqrt(1 - freq.^-2);
figure('Name','d(f/fc) versus f/fc','NumberTitle','off');
plot(freq, dffc);
grid on;
title('d(f/fc)');
xlabel('f/fc');
ylabel('d(f/fc)');

%% Question 2 
%TE and TM modes have same cutoff frequencies for rectangular waveguides
%calculating the cutoff frequencies for air
er = 1;
e0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c = 3e8;
%a and b are the dimensions of the waveguide
a = 0.03;
b = 0.015;
%used to index the modes
m = 0:2;
n = 0:2;

[M,N] = meshgrid(m,n);
%both from formula packet
kc = sqrt((M*pi/a).^2 + (N*pi/b).^2);
%display frequency matrix rather than in a table
fc = kc/(2*pi*sqrt(mu0*e0*er))

%the degenerate modes are TE01 and TM01, and TE10 with itself(?) as they have the same cutoff
%frequency

%for the cutoff frequency for the dominant mode, which is TE10, in AIR, determine which of
%these modes are evanescent
f = 2.5*fc(1,2);
%We can look at the matrix of frequencies and compare each frequency to this one to see which modes 
%are the evanescent modes: TE21, TM21, TE12, TM12, TE22, TM 22, TE 02, not
%in a table

%in AIR find Zte or Ztm as appropriate even for the evanescent case
%from the formulas in the formula packet, we arrive at
eta = sqrt(mu0/(e0*er));
%displaying the two impedance matrices rather than in a table
Zte = eta./sqrt(1-(fc/f).^2) 
Ztm = eta*sqrt(1-(fc/f).^2)

%in AIR for the non evanescent modes, find lambdag,phase velocity, and
%group velocity
v = 1/sqrt(mu0*e0*er);
k = (2*pi*f)/v;
%from formula packet
beta = real(k*sqrt(1-(fc/f).^2))';
beta = beta(1:5);
%from formula packet, displaying them rather than putting them in a table
lambdag = 2*pi./beta 
phasevelocity = v./sqrt(1-(fc/f).^2)
groupvelocity = v*sqrt(1-(fc/f).^2)

%in AIR for each of the evanescent modes, find alpha in dB/m
%from the formula sheet, we arrive at
alpha = imag(kc.*sqrt(1-(fc/f).^2));
%picking out the evanescent modes, displayed rather than in table
alpha = alpha(alpha ~= 0)
%% Still question 2 but sectioned off so you can see the values for air in the first section and dielectric here

%in DIELECTRIC calculate the cutoff frequencies
er = 4;

%same a, b, m0, e0,m, and n as from before, so the only thing we needed to
%change in fc was adding er
%display frequency matrix rather than in a table
fc = kc/(2*pi*sqrt(mu0*e0*er))
%the degenerate modes are TE01 and TM01, and TE10 with itself(?) as they have the same cutoff
%frequency

%for the cutoff frequency for the dominant mode, which is TE10, in DIELECTRIC, determine which of
%these modes are evanescent 
f = 2.5*fc(1,2);
%We can look at the matrix of frequencies and compare each frequency to this one to see which modes 
%are the evanescent modes: TE21, TM21, TE12, TM12, TE22, TM 22, TE 02, not
%in a table

%in DIELECTRIC find Zte or Ztm as appropriate even for the evanescent case
%from the formulas in the formula packet, we arrive at
eta = sqrt(mu0/(e0*er));
%displaying the two impedance matrices rather than in a table
Zte = eta./sqrt(1-(fc/f).^2) 
Ztm = eta*sqrt(1-(fc/f).^2)

%in DIELECTRIC for the non evanescent modes, find lambdag,phase velocity, and
%group velocity
v = 1/sqrt(mu0*e0*er);
k = 2*pi*f/v;
%all from formula packet
beta = real(k*sqrt(1-(fc/f).^2))'
beta = beta(1:5);
%displaying them rather than putting them in a table
lambdag = 2*pi./beta 
phasevelocity = v./sqrt(1-(fc/f).^2)
groupvelocity = v*sqrt(1-(fc/f).^2)

%in DIELECTRIC for each of the evanescent modes, find alpha in dB/m
%from the formula sheet, we arrive at
alpha = imag(kc.*sqrt(1-(fc/f).^2));
%picking out the evanescent modes, displayed rather than in table
alpha = alpha(alpha ~= 0)


%% Question 3
%radius for circular waveguide
a = 0.01;
%positive roots of the derivative of the nth order bessel function for
%TEnm, scaling by 1/a to get a matrix of kc
TEkc = [3.832 7.016 10.174;1.841 5.331 8.536; 3.054 6.706 9.970]/a;
%positive roots of the nth order bessel function of TMnm, scaling by 1/a to
%get a matrix of kc
TMkc = [2.405 5.520 8.650; 3.832 7.016 10.174; 5.135 8.417 11.620]/a;

TEfc = TEkc/(2*pi*sqrt(mu0*e0*er))
TMfc = TMkc/(2*pi*sqrt(mu0*e0*er))
%assuming an operating frequency 10 percent above the largest of these
%cutoffs, at this one frequency, list the guide wavelength, wave impedance,
%and group velocity for the modes
f = 10.0562e9;
fc = [TEfc; TMfc];
%putting them in order for ease
fc = sort(fc(1:5,1));
v = 1/sqrt(mu0*e0*er);
k = (2*pi*f)/v;
beta = k.*sqrt(1-(fc/f).^2);
lambdag = 2*pi./beta;
phasevelocity = v./sqrt(1-(fc/f).^2);
groupvelocity = v*sqrt(1-(fc/f).^2);
%dominant mode is TE11, degenerate modes are TE01 and TM11, see paper for
%tabulated frequencies and such.
