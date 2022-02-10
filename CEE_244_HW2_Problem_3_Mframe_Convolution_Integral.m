%% CEE 244 Homework 2 - Problem 3 Convolution Integrals
%
% Ray Abbiatici
% Version 1.0/RJA/25-Jan-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

clear all
clc

%% Define System Properties

EIcol = 2*10^8;              %   [kN-mm]
L1 = 5*1000;                 %   [mm] 
L2 = 2.5*1000;               %   [mm]
k1 = 12*EIcol/L1^3;          %   [kN/mm]
k2 = 12*EIcol/L2^3;          %   [kN/mm]

%% Compute the PDF of the Moment Demand 

syms x
mup = [0:1:100];

for i = 1:length(mup)

    %% Compute Moment Demands 
    Mu1 = k1/(k1+k2) * L1/2 * mup(i);   %   [kN-mm]
    Mu2 = k2/(k1+k2) * L2/2 * mup(i);   %   [kN-mm]

    a_m1 = 4*pi()/(Mu1*sqrt(6));
    a_m2 = 4*pi()/(Mu2*sqrt(6));

    u_m1 = Mu1*(1-0.577*sqrt(6)/(4*pi));
    u_m2 = Mu2*(1-0.577*sqrt(6)/(4*pi));
    
    fMu_M1 = a_m1*exp(-a_m1*(x-u_m1)-exp(-a_m1*(x-u_m1)));
    fMu_M2 = a_m2*exp(-a_m2*(x-u_m2)-exp(-a_m2*(x-u_m2)));
    
    FMy = (x-5000)/10000;
    pf1(i) = int(fMu_M1*FMy,x,5000,15000)+int(fMu_M1*1,x,15000,inf);
    pf2(i) = int(fMu_M2*FMy,x,5000,15000)+int(fMu_M2*1,x,15000,inf);
end

plot(mup,pf1,mup,pf2)
yline(0.05)
xlabel('Mean of P [kips]')
ylabel('Probability of Failure,pf')
legend('Column A','Column B')

idy1 = find(pf1 == 0.05);
mup05_1 = mup(idy1);
u05_1 = k1/(k1+k2) * L1/2 * mup(idy1)*(1-0.577*sqrt(6)/(4*pi))

idy2 = find(pf2 == 0.05);
mup05_2 = mup(idy2);
u05_2 = k2/(k1+k2) * L2/2 * mup(idy2)*(1-0.577*sqrt(6)/(4*pi))