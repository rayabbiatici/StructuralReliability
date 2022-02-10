% CEE 246 Homework 3 - Reliability Indices - P1&2
%
% Ray Abbiatici
% Version 1.1/RJA/3-Feb-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

clear
clc

%% Define The RV Matricies

mPd = 1500; % [k]
mVe = 150;  % [k]
mfy = 50;   % [ksi]

M = [mPd; mVe; mfy];

%% Define the Coefficient of Variation Matrix

dPd = 0.1;
dVe = 0.25;
dfy = 0.1;

d = [dPd; dVe; dfy];

%% Compute the Standard Deviations

sd = M .* d;

%% Correlation Coefficient Matrix
%    1 - P, 2 - V, 3 - fy

rho = [1.0 0.75 0.0; 0.75 1.0 0.0; 0.0 0.0 1.0];

%% Compute Moment Demands, Flexural Capacity and Yield Force
syms Pd Ve fy
RV = [Pd; Ve; fy];

E = 29000;          %  [ksi]
L = 120;            %  [in]
b = 10;             %  [in]
h = 20;             %  [in]
I = 1/12*b*h^3;     %  [in^4]

Md = Pd * Ve * L^3 / (3 * E * I);
Me = Ve * L;
My = b * h^2 / 6 * fy;
Py = b * h * fy;

%% Define the Limit State Functions

g1 = 1 - (Pd/Py)^2 - ((Md + Me)/My)^2;
g2 = My^2 * Py^2 - Pd^2 * My^2 - Py^2 * (Md + Me)^2;

%% Compute the Reliability index for g1

root = 0;
g1_m = double(subs(g1,{Pd,Ve,fy},{M(1),M(2),M(3)}));

for i = 1:length(M)
    for j = 1:length(M)
        
        corr = rho(i,j);
        sigma_i = sd(i);
        sigma_j = sd(j);
               
        dg1_dxi = diff(g1,RV(i));
        dg1_dxj = diff(g1,RV(j));

        dg1_dxi_m = double(subs(dg1_dxi,{Pd,Ve,fy},{M(1),M(2),M(3)}));
        dg1_dxj_m = double(subs(dg1_dxj,{Pd,Ve,fy},{M(1),M(2),M(3)}));
        
        d_root = dg1_dxi_m * dg1_dxj_m * corr * sigma_i * sigma_j;
        root = root + d_root;
    end
end

beta1 = g1_m/sqrt(root);

%% Compute the Reliability index for g2

root = 0;
g2_m = double(subs(g2,{Pd,Ve,fy},{M(1),M(2),M(3)}));

for i = 1:length(M)
    for j = 1:length(M)
        
        corr = rho(i,j);
        sigma_i = sd(i);
        sigma_j = sd(j);
               
        dg2_dxi = diff(g2,RV(i));
        dg2_dxj = diff(g2,RV(j));

        dg2_dxi_m = double(subs(dg2_dxi,{Pd,Ve,fy},{M(1),M(2),M(3)}));
        dg2_dxj_m = double(subs(dg2_dxj,{Pd,Ve,fy},{M(1),M(2),M(3)}));

        d_root = dg2_dxi_m * dg2_dxj_m * corr * sigma_i * sigma_j;
        root = root + d_root;
    end
end

beta2 = g2_m/sqrt(root);

%% Begin Problem 2

%% Define the LSF and key matrices
RV = [Pd; Ve; fy];
M = [mPd; mVe; mfy];

E = 29000;          %  [ksi]
L = 120;            %  [in]
b = 10;             %  [in]
h = 20;             %  [in]
I = 1/12*b*h^3;     %  [in^4]

Md = Pd * Ve * L^3 / (3 * E * I);
Me = Ve * L;
My = b * h^2 / 6 * fy;
Py = b * h * fy;

g = 1 - (Pd/Py)^2 - ((Md + Me)/My)^2;

D = zeros(length(sd));

for i = 1:length(sd)
    D(i,i) = sd(i);
end

L = (chol(rho))';

%% Compute the Gradient of h(u) - the LSF in U Space

for i = 1:length(RV)
    grad_g(i) = diff(g,RV(i));
end


%% Begin HL Iteration
DL = D*L;

e1 = 1;
e2 = 1;

k = 1;

x(:,k) = [1500;150;50];
u(:,k) = inv(L)*inv(D)*(x(:,k)-M);

grad_g_i = double(subs(grad_g,{Pd,Ve,fy},{x(1,k),x(2,k),x(3,k)}));
h = double(subs(g,{Pd,Ve,fy},{x(1,k),x(2,k),x(3,k)}));
grad_h = (DL)' * grad_g_i';

alpha(:,k) = -grad_h./norm(grad_h);

beta(:,k) = alpha(:,k)' * u(:,k);

u(:,k+1) = alpha(:,k) * (beta(:,k) + h/norm(grad_h));

x(:,k+1) = D * L * u(:,k+1) + M;

while (e1 > 10^-3 && e2 > 10^ -3)
    
    k = k + 1;
    
    grad_g_i = double(subs(grad_g,{Pd,Ve,fy},{x(1,k),x(2,k),x(3,k)}));
        h = double(subs(g,{Pd,Ve,fy},{x(1,k),x(2,k),x(3,k)}));
    
    grad_h = (D*L)' * grad_g_i';
    
    alpha(:,k) = -grad_h ./ norm(grad_h);
    
    beta(:,k) = alpha(:,k)' * u(:,k);
    
    u(:,k+1) = alpha(:,k) * (beta(:,k) + h/norm(grad_h));
    
    x(:,k+1) = D * L * u(:,k+1) + M;
    
    e1 = abs(beta(:,k) - beta(:,k-1));
    e2 = abs(h);
    
end

% Compute Importance Vector
gamma = inv(L)'*alpha(:,k-1)/norm(inv(L)'*alpha(:,k-1));

