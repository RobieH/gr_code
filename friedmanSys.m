function [A, T, R, W] = friedmanSys(Or,Om,Ov, Trange, eta, maxcount)
H = 13.9;
G = 7.426E-28;
tol = 1.05;
%rhoC =  8*pi*G/(3*H^2);
rhoC = 1;
Omegas = [Or, Om, Ov];

%solve the system of coupled DEs using 4th order Runge-Kutta in forward and
%reverse

Af = [1];
Rf = [rhoC];
Tf = [0];
Wf = [w(1,Omegas)];

i = 2;
count = 0;
step = eta;
deltaTf = 0;
while deltaTf < Trange
    %calculate RK parameters for the DEs.  Use k for rho de and m for a de
    k1 = rhohelper(Af(i-1), Rf(i-1),Omegas);
    m1 = ahelper(Af(i-1), Rf(i-1));
    
    k2 = rhohelper(Af(i-1)+step*m1/2, Rf(i-1)+step*k1/2, Omegas);
    m2 = ahelper(Af(i-1)+step*m1/2, Rf(i-1)+step*k1/2);
    
    k3 = rhohelper(Af(i-1)+step*m2/2, Rf(i-1)+step*k2/2, Omegas);
    m3 = ahelper(Af(i-1)+step*m2/2, Rf(i-1)+step*k2/2);
    
    k4 = rhohelper(Af(i-1)+step*m3, Rf(i-1)+step*k3, Omegas);
    m4 = ahelper(Af(i-1)+step*m3, Rf(i-1)+step*k3);
    
    Af(i) = Af(i-1) + step/6*(m1+2*m2+2*m3+m4);
    Rf(i) = Rf(i-1) + step/6*(k1+2*k2+2*k3+k4);
    Tf(i) = Tf(i-1) + step*Af(i)*H;
    Wf(i) = w(Af(i), Omegas);
    
    diffA = Af(i)/Af(i-1);
    diffR = Rf(i)/Rf(i-1);
    
    if (diffA > tol || diffR > tol) && count < maxcount
        step = step/10;
        count = count + 1;
    elseif (diffA > tol || diffR > tol) && count >= maxcount
        break;
    else 
        deltaTf = Tf(i);
        i = i + 1;
    end
        
    

    
    
end

% now backward Runge-Kutta
Ab = [1];
Rb = [rhoC];
Tb = [0];
Wb = [];
step = -eta; %iterate backwards
i = 2;
deltaTb = 0;
count = 0;

while deltaTb < Trange

    %calculate RK parameters for the DEs.  Use k for rho de and m for a de
    k1 = rhohelper(Ab(i-1), Rb(i-1),Omegas);
    m1 = ahelper(Ab(i-1), Rb(i-1));
    
    k2 = rhohelper(Ab(i-1)+step*m1/2, Rb(i-1)+step*k1/2, Omegas);
    m2 = ahelper(Ab(i-1)+step*m1/2, Rb(i-1)+step*k1/2);
    
    k3 = rhohelper(Ab(i-1)+step*m2/2, Rb(i-1)+step*k2/2, Omegas);
    m3 = ahelper(Ab(i-1)+step*m2/2, Rb(i-1)+step*k2/2);
    
    k4 = rhohelper(Ab(i-1)+step*m3, Rb(i-1)+step*k3, Omegas);
    m4 = ahelper(Ab(i-1)+step*m3, Rb(i-1)+step*k3);
    
    Ab(i) = Ab(i-1) + step/6*(m1+2*m2+2*m3+m4);
    Rb(i) = Rb(i-1) + step/6*(k1+2*k2+2*k3+k4);
    Tb(i) = Tb(i-1) + step*Ab(i)*H;
    Wb(i) = w(Ab(i), Omegas);
    
    diffA = Ab(i)/Ab(i-1);
    diffR = Rb(i)/Rb(i-1);
    
    if (diffA > tol || diffR > tol) && count < maxcount
        step = step/10;
        count = count + 1;
    elseif (diffA > tol || diffR > tol) && count >= maxcount
        break;
    else
        deltaTb = -Tb(i); 
        i = i + 1;
 
    end 
end
    
A = [Ab(end:-1:1), Af(2:length(Af))]; %reverse vector and splice together
R = [Rb(end:-1:1), Rf(2:length(Rf))];
T = [Tb(end:-1:1), Tf(2:length(Tf))];
W = [Wb(end:-1:2), Wf(1:length(Wf))];

%there may be NaNs at the beginning, due to a rapidly going to zero.  here
%we check for them and delete them if present.

removeA = isnan(A) + isinf(A);
removeR = isnan(R) + isinf(R);
removeT = isnan(T) + isinf(T);
removeW = isnan(W) + isinf(W);

remove = removeA + removeR + removeT + removeW;
toDelete = [remove>0];

A(toDelete == 1) = [];
R(toDelete == 1) = [];
T(toDelete == 1) = [];
W(toDelete == 1) = [];


plot(T,A, T, R);    

    