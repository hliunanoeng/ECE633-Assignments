%graphene band structure using the tight-binding method

a_cc = 1.42e-10; %[m]
eps = 0; %[eV]
t = -3; %[eV]
a = (3/2)*a_cc; %[m]
b = ((3^(1/2))/2)*a_cc; %[m}

a1 = [a;b];
a2 = [a;-b];

l = 101;
k_x = linspace(-pi/a, pi/a, l);
k_y = linspace(-pi/b, pi/b, l);

[k_X,k_Y] = meshgrid(k_x,k_y);

E1 = [];
E2 = [];

for ii=1:l
    for jj = 1:l
        h = [eps,0;0,eps];
        k = [k_x(ii);k_y(jj)];
        h_0 = t*(1+exp(i*dot(k,a1))+exp(i*dot(k,a2)));
        h(1,2) =  h_0';
        h(2,1) = h_0;
        [V,D] = eig(h);
        d = diag(D);
        E1(ii,jj) = d(1);
        E2(ii,jj) = d(2);
    end
end

surf(k_X,k_Y,E1);
hold on;
surf(k_X,k_Y,E2);
xlabel('k_{x} [1/m]');
ylabel('k_{y} [1/m]');
zlabel('E [eV]');
