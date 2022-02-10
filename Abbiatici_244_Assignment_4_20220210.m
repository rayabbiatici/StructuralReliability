% CEE 244 Assignment 4 - Monte Carlo Methods
%
% Ray Abbiatici
% Version 1.0/RJA/10-Feb-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%
%
%% Dictionary of Variables
%
% Inputs:
%   lw = Length of Wall [in]
%   hoops = Confining Reinf. Indicator [0 for overlap, 1 for single]
%   tw = Wall thickness [in]
%   a = Empirical Constant for Spectral Acceleration Computation
%   b = Empirical Constrant for Spectral Acceleration Computation
%   n = Number of samples to be realized
%   CV_x = Coefficient of Variation of RV x
%
% Intermediate Variables:
%   c = Neutral Axis Depth [in]
%   bw = Wall dimension parallel to bending axis [in]
%   DBE = First mode spectral acceleration for DBE [in/sec^2]
%   MCE = First mode spectral acceleration for DBE [ in/sec^2]
%   mean_XXX = Mean value of random variable XXX [varies]
%   std_XXX = Standard Deviation of random variable XXX [varies]
%   Y = Vector of realizations of CDF of random variable in question
%   XX_bar = Vector of realizations of random variable XX [varies]
%   vmax = Maximum Shear Stress [kips]
%
% Outputs:
%   pf = Probability of Failure due to MCE Shaking
%   pf_cond = Probability of Failure condition on DBE shaking exceeding
%       1/2 of drift capacity
%   ans_a = Vector of estimates of pf
%   ans_b = Vector of estimates of pf_con

for j = 1:10

%% Define Inputs & Intermediate Parameters

lw = 14*12; % [in]
hoops = 1;  % [0 for overlap, 1 for single]
tw = 14;    % [in]
bw = tw;    % [in]

%% Sample the CDF of each Random Variable
n = 1000000;            % [Samples realized]
Y_fc = rand(n,1);
Y_vmax = rand(n,1);
Y_c = rand(n,1);
Y_theta_c = rand(n,1);
Y_theta_max_DBE = rand(n,1);
Y_theta_max_MCE = rand(n,1);

%% Provide Realizations of distributed parameters

% f'c - Concrete Compressive Strength (Exponential)
mean_fc = 1/0.2; % [ksi]
std_fc = 1/0.2;  % [ksi]

fc_bar = expinv(Y_fc,mean_fc); % [ksi]

% vmax - Maximum Shear Stress (Rayleigh)
mean_vmax = 6*sqrt(fc_bar*1000);
alpha = mean_vmax*2/sqrt(pi);
        
B = alpha/sqrt(2);
vmax_bar = raylinv(Y_vmax,B);

% c - Depth of Neutral Axis 
mean_c = 0.25*lw; % [in]
CV_c = 0.2;

std_ln_c = sqrt(log(CV_c^2+1));
mean_ln_c = log(mean_c) - 1/2*(std_ln_c)^2;

c_bar = logninv(Y_c,mean_ln_c,std_ln_c);

% theta_c - Drift Capacity (Normal)
lambda_b = lw*c_bar/bw^2;

if hoops == 1
    alpha_hoop = 45;
else
    alpha_hoop = 60;
end

mean_theta_c = 3.85 - lambda_b./alpha_hoop - vmax_bar./(10*sqrt(fc_bar*1000));
CV_theta_c = 0.15;
std_theta_c = mean_theta_c * CV_theta_c;

theta_c_bar = norminv(Y_theta_c,mean_theta_c,std_theta_c);

% theta_max - Drift Demand (Lognormal)
a = 0.9;
b = 1.1;

SaT1_DBE = 1;
SaT1_MCE = 1.5;

med_theta_max_DBE = a*SaT1_DBE^b;
med_theta_max_MCE = a*SaT1_MCE^b;

mean_ln_theta_max_DBE = log(med_theta_max_DBE);
mean_ln_theta_max_MCE = log(med_theta_max_MCE);

std_ln_theta_max = 0.4;

theta_max_DBE_bar = logninv(Y_theta_max_DBE,mean_ln_theta_max_DBE, ...
    std_ln_theta_max);

theta_max_MCE_bar = logninv(Y_theta_max_MCE,mean_ln_theta_max_MCE, ...
    std_ln_theta_max);

%% Develop the Limit State Function

g = theta_c_bar - theta_max_MCE_bar;
count = 0;

%% Compute the probability of failure

for i = 1:length(g)
    
    if g(i) <= 0
        count = count + 1;
    end
end

pf = count/length(g);
ans_a(j) = pf;


%% Compute the probability of failure provided that 1.2 of the capacity
%   is exceeded

count = 0;
space = 0;
for i = 1:length(g)
    if theta_max_DBE_bar(i) > 0.5*theta_c_bar(i) && g(i) <= 0
        space = space+1;
        count = count + 1;
    

    elseif theta_max_DBE_bar(i) > 0.5*theta_c_bar(i) && g(i) > 0
        space = space + 1; 
    end
end

pf_cond = count/space;
ans_b(j) = pf_cond;

end 

%% Select bw to limit pf to 5%
lw = 14*12; % [in]
hoops = 1;  % [0 for overlap, 1 for single]
tw = 14;    % [in]
bw = tw;  % [in]

while pf > 0.05

%% Sample the CDF of each Random Variable
n = 1000000;            % [Samples realized]
Y_fc = rand(n,1);
Y_vmax = rand(n,1);
Y_c = rand(n,1);
Y_theta_c = rand(n,1);
Y_theta_max_DBE = rand(n,1);
Y_theta_max_MCE = rand(n,1);

%% Provide Realizations of distributed parameters

% f'c - Concrete Compressive Strength (Exponential)
mean_fc = 1/0.2; % [ksi]
std_fc = 1/0.2;  % [ksi]

fc_bar = expinv(Y_fc,mean_fc); % [ksi]

% vmax - Maximum Shear Stress (Rayleigh)
mean_vmax = 6*sqrt(fc_bar*1000);
alpha = mean_vmax*2/sqrt(pi);
        
B = alpha/sqrt(2);
vmax_bar = raylinv(Y_vmax,B);

% c - Depth of Neutral Axis 
mean_c = 0.25*lw; % [in]
CV_c = 0.2;

std_ln_c = sqrt(log(CV_c^2+1));
mean_ln_c = log(mean_c) - 1/2*(std_ln_c)^2;

c_bar = logninv(Y_c,mean_ln_c,std_ln_c);

% theta_c - Drift Capacity (Normal)
lambda_b = lw*c_bar/bw^2;

if hoops == 1
    alpha_hoop = 45;
else
    alpha_hoop = 60;
end

mean_theta_c = 3.85 - lambda_b./alpha_hoop - vmax_bar./(10*sqrt(fc_bar*1000));
CV_theta_c = 0.15;
std_theta_c = mean_theta_c * CV_theta_c;

theta_c_bar = norminv(Y_theta_c,mean_theta_c,std_theta_c);

% theta_max - Drift Demand (Lognormal)
a = 0.9;
b = 1.1;

SaT1_DBE = 1;
SaT1_MCE = 1.5;

med_theta_max_DBE = a*SaT1_DBE^b;
med_theta_max_MCE = a*SaT1_MCE^b;

mean_ln_theta_max_DBE = log(med_theta_max_DBE);
mean_ln_theta_max_MCE = log(med_theta_max_MCE);

std_ln_theta_max = 0.4;

theta_max_DBE_bar = logninv(Y_theta_max_DBE,mean_ln_theta_max_DBE, ...
    std_ln_theta_max);

theta_max_MCE_bar = logninv(Y_theta_max_MCE,mean_ln_theta_max_MCE, ...
    std_ln_theta_max);

%% Develop the Limit State Function

g = theta_c_bar - theta_max_MCE_bar;
count = 0;

%% Compute the probability of failure

for i = 1:length(g)
    
    if g(i) <= 0
        count = count + 1;
    end
end

pf = count/length(g);
bw = bw + 0.25;

end

disp(bw-0.25)


