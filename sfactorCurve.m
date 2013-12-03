function [A, T] = sfactor(Or,Om,Ov, Trange, eta, maxCount)
%   a(t) Creates a plot of a(t) vs. t (in Gy) based on Omega parameters.
%   Implements a fourth-order Runge-Kutta integrator with the ability to 
%   adapt the step size in a simple way for better accuracy.
%
%   The output is the scale factor of the spatial part of the universe, 
%   showing the size of the universe relative to its size today (i.e.
%   a(now) = 1).  
%
%   Requires a helperCurve function (denoted helperCurve below).  This function
%   should be a function of the omegas and a.  It can be chosen based upon
%   suitable choice of the equation of state parameter, w, in the solution
%   of the Friedmann equation.
%
%   INPUTS
%   Or - radiation density (0.000084 for our universe)
%   Om - mass denisity (0.272 for our universe)
%   Ov - Vacuum energy density (0.728 for our universe)
%   Trange - range of times for the a(t) plot.  
%   eta - initial step size for iterations (0.01 seems to be good)
%   maxCount - the maximum number of times the program should fail to
%       produce a reasonable output before it is terminated. 
%
%   OUTPUTS
%   A - A vector containing the computed a(t) values.
%   T - A vector containing the computed t values.


Omegas = [Or, Om, Ov]; %vector of omegas for helperCurve function

%first calculate points in time forward from today

Af = [1];
Tf = [0];
i = 2;
count = 0;
step = eta;
deltaTf = 0;
while deltaTf < Trange
    %calculate RK parameters
    k1 = helperCurve(Af(i-1), Omegas);
    k2 = helperCurve(Af(i-1)+step/2*k1, Omegas);
    k3 = helperCurve(Af(i-1)+step/2*k2, Omegas);
    k4 = helperCurve(Af(i-1)+step*k3, Omegas);
    
    if ~isreal([k1,k2,k3,k4])
        break
    end
    
    %evaluate new A, T values
    Af(i) = Af(i-1) + step/6*(k1+2*k2+2*k3+k4);
    Tf(i) = Tf(i-1)+Af(i)*step*13.9;
    
    %check to make sure step size is still appropriate
    diff = Af(i)/Af(i-1);
    deltaTf = Tf(i);
    %determine if the value of A(i) is ok
    if diff > 1.1 && count < maxCount %1.1 was determined by trial/error
        step = step/10
        if step < 1000*eps, break, end %don't want step size < eps
        count = count + 1
    elseif diff > 1.1 && count >= maxCount
        break
    else 
        i = i+1;
    end
end

%Now calculate the evolution of the early universe
Ab = [1];
Tb = [0];

step = -eta; %iterate backwards
i = 2;
deltaTb = 0;

while deltaTb < Trange
    
    k1 = helperCurve(Ab(i-1), Omegas);
    k2 = helperCurve(Ab(i-1)+step/2*k1, Omegas);
    k3 = helperCurve(Ab(i-1)+step/2*k2, Omegas);
    k4 = helperCurve(Ab(i-1)+step*k3, Omegas);
    
    if ~isreal([k1,k2,k3,k4])
        break
    end
    
    Ab(i) = Ab(i-1) + step/6*(k1+2*k2+2*k3+k4);
    Tb(i) = Tb(i-1)+Ab(i)*step*13.9;

    diff = Ab(i)/Ab(i-1);
    deltaTb = -Tb(i);
    if diff > 1.1 && count < maxCount
        step = step/10;
        if step < 1000*eps, break, end
        count = count + 1;
    elseif diff > 1.1 && count >= maxCount
        break
    elseif i > 400
        break
    else 
        i = i+1;
    end
end

A = [Ab(end:-1:1), Af(2:length(Af))]; %reverse vector and splice together
T = [Tb(end:-1:1), Tf(2:length(Tf))];
plot(T,A);
end




