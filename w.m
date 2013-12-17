function w = w(a, Omega)

% k1 = 0.1;
% k2 =5;
% 
% val1 = - (1+exp(-k2*(0.728-0.272*a^(-3))))^-1;
% val2 = - 2/3*(1+exp(k1*(0.272*a^(-3)-0.000085*a^(-4))))^-1;
% w = 1/3 + val1 + val2;

ar_max = Omega(1)/Omega(2);

%matter epoch 

am_max = (Omega(2)/Omega(3))^(1/3);

%determine what value of w to use (the equivalent of a Heaviside function)

if a < ar_max
    w = 1/3;
elseif a < am_max && a > ar_max
    w = 0;
% elseif a >= 1
%    w = -1.19;
else
    w = -1; 
end



