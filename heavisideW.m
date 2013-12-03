function val = heavisideW(a, Omega)

%helper function for sfactor

%first we need to determine the epochs

%radiation epoch

ar_max = Omega(1)/Omega(2);

%matter epoch 

am_max = (Omega(2)/Omega(3))^(1/3);

%determine what value of w to use (the equivalent of a Heaviside function)

if a < ar_max
    w = 1/3;
    Om = Omega(1);
elseif a < am_max && a > ar_max
    w = 0;
    Om = Omega(2);
elseif a >= 1
    w = -1.19;
    Om = 1;
else
    w = -1; 
    Om = Omega(3);
end


val = sqrt(Om)*a.^((-3*w+1)/2);