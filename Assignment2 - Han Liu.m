clear all;  % to clear all parameters, if any
% close all;  % to close all figures, if any

dE = 1e-3; % integration step [eV]
gamma_1 = 0.2; %[eV] 
gamma_2 = 0.2; %[eV] 
mu_equilibrium = 0; %[eV]

%populate the arrays
l_V = [0:0.001:0.5]; %[V]
l_I = [];
for n = 1:length(l_V)
    v = l_V(n);
    i = current(v,mu_equilibrium,gamma_1,gamma_2,dE);
    l_I(n) = i;
end

%plot the data
plot(l_V,l_I.*1.0e6,'b-','linewidth',2)
xlabel('Voltage (V)');
ylabel('Current (\muA)');
set(gca,'xtick',[0,0.1,0.2,0.3,0.4,0.5]);

function [df,E] = dFermi(mu1,mu2,dE)
    % mu1 > mu2 [eV]
    E = [mu2-1:dE:mu1+1]; % [eV]
    q = 1.6e-19; %[C]
    kB = 1.38e-23;      % [J/K]
    kB = kB/q;          % [eV/K]
    T = 300;            % [K]
    f1 = 1 ./ (1+exp((E-mu1)./(kB*T)));
    f2 = 1 ./ (1+exp((E-mu2)./(kB*T)));
    df = f1 - f2;
end

function D = DOS(E)
    E = E>= -0.1; %[eV]
    D = E.* 0.5; %[eV^-1]
end

function I = current(V,mu_eq,gamma1,gamma2,dE)
    q = 1.6e-19; %[C]
    h_bar = 6.626e-34/(2*pi); %[J*s]
    h_bar = h_bar*6.242e18; %[eV*s]
    delta_e = q*V*6.242e18; %[eV]
    mu1 = mu_eq + (delta_e/2); %[eV]
    mu2 = mu_eq - (delta_e/2); %[eV]
    [df,E] = dFermi(mu1,mu2,dE);
    gamma_factor = (gamma1*gamma2)/(gamma1+gamma2); %[eV]
    D = DOS(E); %[eV^-1]
    I = ((2*q)/h_bar).*sum(D.*gamma_factor.*df.*dE); %[A]
end