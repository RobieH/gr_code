function rhoval = rhohelper(a, rho, Omegas)

H = 7.61E-27;
G = 7.426E-28;

%rhoC = 8*pi*G/(3*H^2);
rhoC = 1;
rhoval = a*rhoC*(-3*(1+w(a,Omegas))*rho^(3/2));
