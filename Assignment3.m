clear all;  % to clear all parameters, if any
% close all;  % to close all figures, if any

dE = 1e-3; % integration step [eV]
gamma_1 = 0.2; %[eV] 
gamma_2 = 0.2; %[eV] 
mu_equilibrium = 0; %[eV]
U_S = 0;%[eV]
mu_1 = mu_equilibrium+U_S;%[eV]

%line 1

U_0 = 1; %[eV]

%populate the arrays
l_V = [0:0.01:0.5]; %[V]
l_I_1 = [];

for n = 1:length(l_V)
    v = l_V(n);%[V]
    U_D = (-1)*v;%[eV]
    mu_2 = mu_equilibrium+U_D;%[eV]
    U = SCS(v,U_0,mu_1,mu_2,gamma_1,gamma_2,dE);
    i = current(mu_1,mu_2,gamma_1,gamma_2,U,dE);
    l_I_1(n) = i;
end

%line 2

U_0 = 3; %[eV]

%populate the arrays
l_V = [0:0.01:0.5]; %[V]
l_I_2 = [];

for n = 1:length(l_V)
    v = l_V(n);%[V]
    U_D = (-1)*v;%[eV]
    mu_2 = mu_equilibrium+U_D;%[eV]
    U = SCS(v,U_0,mu_1,mu_2,gamma_1,gamma_2,dE);
    i = current(mu_1,mu_2,gamma_1,gamma_2,U,dE);
    l_I_2(n) = i;
end

%plot the data
plot(l_V,l_I_1.*1.0e6,'b-','linewidth',2)
hold on
plot(l_V,l_I_2.*1.0e6,'r-.','linewidth',2)
hold off
xlabel('Voltage (V)');
ylabel('Current (\muA)');
set(gca,'xtick',[0,0.1,0.2,0.3,0.4,0.5]);
leg = legend('U_{L}=1eV','U_{L}=3eV');
leg.Location = 'northwest';

function [df,E,f1,f2] = dFermi(mu1,mu2,dE)
    % mu1 > mu2 [eV]
    E = [-1:dE:1]; % [eV]
    q = 1.6e-19; %[C]
    kB = 1.38e-23;      % [J/K]
    kB = kB/q;          % [eV/K]
    T = 300;            % [K]
    f1 = 1 ./ (1+exp((E-mu1)./(kB*T)));
    f2 = 1 ./ (1+exp((E-mu2)./(kB*T)));
    df = f1 - f2;
end

function D = DOS(E)
    E = E>= 0; %[eV]
    D = E.* 0.5; %[eV^-1]
end

function I = current(mu1,mu2,gamma1,gamma2,U,dE)
    q = 1.6e-19; %[C]
    h_bar = 6.626e-34/(2*pi); %[J*s]
    h_bar = h_bar*6.242e18; %[eV*s]
    [df,E] = dFermi(mu1,mu2,dE);
    gamma_factor = (gamma1*gamma2)/(gamma1+gamma2); %[eV]
    D = DOS(E-U); %[eV^-1]
    I = ((2*q)/h_bar).*sum(D.*gamma_factor.*df.*dE); %[A]
end

function U = U_sim(U_L,U_0,N,N_eq)
    U = U_L + U_0*(N-N_eq);%[eV]scalar
end

function N = N_sim(U,mu1,mu2,gamma1,gamma2,dE)
    [df,E,f1,f2]=dFermi(mu1,mu2,dE);
    D = DOS(E-U); %[eV^-1]
    gamma_term = ((gamma1*f1)+(gamma2*f2))/(gamma1+gamma2);
    N = 2*sum((D.*gamma_term)*dE);%scalar
end

function N_eq = N_eq_sim(U_L,mu,dE)
    [df,E,f1,f2]=dFermi(mu,mu,dE);
    D = DOS(E-U_L); %[eV^-1]
    N_eq = 2*sum((D.*f1)*dE); %scalar
end

function final_U = SCS(V,U_0,mu1,mu2,gamma1,gamma2,dE)
    mu = 0;%[eV]
    U_L = 0;
    
    N_eq = N_eq_sim(U_L,mu,dE);%scalar
    
    U_S = 0;%[eV]
    U_D = (-1)*V;%[eV]
    U_L = 0.5*(U_S+U_D);%[eV]
    
    U_init = U_sim(U_L,U_0,N_eq,N_eq);%[eV]scalar
    N = N_sim(U_init,mu1,mu2,gamma1,gamma2,dE);%scalar
    U_prev = U_sim(U_L,U_0,N,N_eq);
    
    diff = 1;%while loop guranteed
    
    while diff>1e-3 %[eV]
        N = N_sim(U_prev,mu1,mu2,gamma1,gamma2,dE);%scalar
        U_curr = U_sim(U_L,U_0,N,N_eq);%[eV]scalar
        diff = abs(U_curr - U_prev);
        %U_prev = U_curr;
        U_prev = 0.5 * U_curr + 0.5 * U_prev;
    end
    
    final_U = U_curr;
    
end