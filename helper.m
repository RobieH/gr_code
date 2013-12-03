function val = helper(a, Omega)

%helper function for sfactor

val = sqrt(Omega(1) + Omega(2).*a + Omega(3).*a^4);
