%% CEE 244 Homework 2 - Problem 4 Convolution Integrals
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
mup = [.2:.2:100];
mmy = 10000;
smy = 10000/sqrt(12);
d = 0.25;

for i = 1:length(mup)

    %% Compute Moment Demands 
    Mu1 = k1/(k1+k2) * L1/2 * mup(i);   %   [kN-mm]
    Mu2 = k2/(k1+k2) * L2/2 * mup(i);   %   [kN-mm]
    
    B1 = (mmy-Mu1)/sqrt((d*Mu1)^2+(d*smy)^2);
    B2 = (mmy-Mu2)/sqrt((d*Mu2)^2+(d*smy)^2);

    pf1(i) = normcdf(-B1);
    pf2(i) = normcdf(-B2);

    lambda_0_1(i) = mmy/Mu1;
    lambda_0_2(i) = mmy/Mu2;
    
end

figure
plot(mup,pf1,mup,pf2)
yline(0.05)
xlabel('Mean of P [kips]')
ylabel('Probability of Failure,pf')
legend('Column A','Column B')

figure
plot(mup,lambda_0_1,mup,lambda_0_2)
xlabel('Mean of P [kips]')
ylabel('Central Safety Factor')
legend('Column A','Column B')


