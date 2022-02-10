% CEE 246 Homework 3 - Reliability Indices - P3
%
% Ray Abbiatici
% Version 1.1/RJA/3-Feb-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

clear
clc

%% Begin Problem 3
% Modify rho assuming Pd and Ve are uncorrolated.
syms Pd Ve fy
rho = [1.0 0 0.0; 0 1.0 0.0; 0.0 0.0 1.0];

%% Define the LSF and key matrices
RV = [Pd; Ve; fy];

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


