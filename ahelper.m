function aval = ahelper(a, rho)

H = 7.61E-27;
G = 7.426E-28;

%rhoC = 8*pi*G/(3*H^2);
rhoC = 1;
aval = a^2*sqrt(rhoC*rho);