function val = helperCurve(a, Omega)

%helper function for sfactor


number = Omega(1) + Omega(2)*a + Omega(3).*a^4+(1-sum(Omega).*a^2);
val = sqrt(number);
