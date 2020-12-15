%disclaimer: the following script was modified from Prof.Yoon's FDM.m code

h = 6.626e-34;      % [J-s]
hbar = h/(2*pi);    % [J-s]
m = 9.11e-31;       % [kg]
q = 1.6e-19;        % [C]
L = 5e-9;          % [m]
a = 1e-10;          % [m]
t0 = hbar^2/(2*m*a^2)/q;  % [eV]
mu = 0.5; % [eV]
N = L/a-1;          % BBC

H = 2*t0*diag(ones(N,1)) -t0*diag(ones(N-1,1),1) -t0*diag(ones(N-1,1),-1);

[V,D] = eig(H);
rho = D - mu.*diag(ones(N,1));

for i=1:N
    rho(i,i) = f_0(rho(i,i));
end

Rho = V * rho * ctranspose(V);
Rho = diag(Rho);
Rho = 2*[0;Rho;0]./a;
x = ((0:N+1)*a);

plot(x*(10^9),Rho, 'b-');
set(gca,'fontsize',20);
xlabel('Position [nm]');
ylabel('Density [/m]');


function r = f_0(E)
    q = 1.6e-19;        % [C]
    kB = 1.38e-23; %[J/K]
    kB = kB/q; %[eV/K]
    T = 300; %[K]
    r = 1/(1+exp((E/(kB*T))));
end