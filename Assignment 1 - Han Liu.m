clear all;  % to clear all parameters, if any
% close all;  % to close all figures, if any


q = 1.6e-19;
E = [-1:0.01:0.5]; % [eV]
mu1 = 0;             % [eV]
mu2 = -0.5;         % [eV]
kB = 1.38e-23;      % [J/K]
kB = kB/q;          % [eV/K]
T = 300;            % [K]

f1 = 1 ./ (1+exp((E-mu1)/(kB*T)));
f2 = 1 ./ (1+exp((E-mu2)/(kB*T)));
f = f1 - f2;

plot(f,E,'b-','linewidth',2); hold on;
set(gca,'fontsize',20);

xlabel('f_{1} - f_{2} (E)');
ylabel('E (eV)');
set(gca,'xtick',[0,0.2,0.4,0.6,0.8,1.0]);
set(gca,'ytick',[-1,-0.5,0,0.5]);



