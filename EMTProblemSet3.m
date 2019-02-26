%% Michael Lendino EMT PSET 3
clc;
clear all;
%% Problem 1
%Write code that takes the source parameters Vg, Zg, the load impedance Zl,
%the characteristic impedance Zo, the normalized length l/lambda as inputs,
%and compute Vin, Iin, Vl, Il, and V(z), I(z), the number of points N for z
%should be an input parameter, and the code should produce plots of the
%magnitude of |V(z)| and |I(z)| as subplots, and the code should compute
%the power input to the transmission line, and the power delivered to the
%load, and check that they match

Vg = 1;
Zg = 3+j*4;
%uncomment which Z1 you want depending on the case you want for the load
%impedance
Zl = 10+j*5;
% Zl = 51-j*1;
Zo = 50;
%Thank you aziza and armaan for helping me decipher/figure out how to
%handle the normalized length and put them in my equations, making sure z=0
%at the load
loverlambda = 1.6;
N = 1000;
dist = linspace(-1,0,N);

%From page 9 of Tx line phasors, computation of Vin and Iin
Zin = Zo*(Zl + j*Zo*tan(2*pi*loverlambda))/(Zo + j*Zl*tan(2*pi*loverlambda));
Vin = Zin/(Zin + Zg)*Vg; 
Iin = Vg/(Zin + Zg);
%page 3
gammaload = (Zl - Zo)/(Zl + Zo);
%page 5
gammaIN = gammaload*exp(-j*4*pi*loverlambda);

Vplus = Vin/(exp(j*2*pi*1.6)*(1+gammaIN));
%Computation of V(z) and I(z)
Vz = Vplus*exp(-j*2*pi*dist).*(1+gammaload*exp(j*4*pi*dist));
Iz = (1/Zo)*(Vplus*exp(-j*2*pi*dist).*(1-gammaload*exp(j*4*pi*dist)));
%Computation of Vload and Iload
Vload = Vz(end);
Iload = Iz(end);
%Plots for |V(z)| and |I(z)|
figure('Name','Magnitude Plots for V and I in the Transmission Line','NumberTitle','off');
subplot(2,1,1);
plot(dist, abs(Vz));
xlabel('z')
ylabel('Volts')
grid on;
title('|V(z)| vs z');
subplot(2,1,2);
plot(dist, abs(Iz));
xlabel('z')
ylabel('Amps')
grid on;
title('|I(z)| vs z');
%Computation of the power input to the transmission line and the power
%delivered to the load
Pin = 0.5*real(Vz(1)*conj(Iz(1)));
Pload = 0.5*real(Vz(end)*conj(Iz(end))); 
%They are the same in both cases because we assumed all transmission lines
%are lossless for this problem, so all of the power is delivered to the
%load