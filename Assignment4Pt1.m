%disclaimer: the following script was modified from Prof.Yoon's FDM.m code

h = 6.626e-34;      % [J-s]
hbar = h/(2*pi);    % [J-s]
m = 9.11e-31;       % [kg]
q = 1.6e-19;        % [C]
L = 5e-9;          % [m]
a = 1e-10;          % [m]
t0 = hbar^2/(2*m*a^2)/q;  % [eV]

N = L/a-1;          % BBC
H = 2*t0*diag(ones(N,1)) -t0*diag(ones(N-1,1),1) -t0*diag(ones(N-1,1),-1);

[V,D] = eig(H);
D = diag(D);

n = [1:25];
E = (((hbar^2)*(pi^2)/(2*m*(L^2)))*(n.^2))/q;

plot(D, 'rx'); hold on;
plot(n,E, 'b-'); hold on;
set(gca,'fontsize',20);
xlim([0,25]);
xlabel('Eigen Number');
ylabel('Eigenvalues [eV]');