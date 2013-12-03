function val = analyticW(a, Omega)

k1 = 10;
k2 = 10;

w = 1/3 - 1/3*(1+exp(-k1*(Omega(2)*a^(-3)-Omega(1)*a^(-4))))^-1 - (1+exp(-k2*(Omega(3)-Omega(2)*a^(-3))))^-1;
Om = Omega(1) - (Omega(1)-Omega(2))*(1+exp(-k1*(Omega(2)*a^(-3)-Omega(1)*a^(-4))))^-1 - (Omega(2)-Omega(3))*(1+exp(-k2*(Omega(3)-Omega(2)*a^(-3))))^-1;
disp([Om,a])
val = sqrt(Om)*a.^((-3*w+1)/2);