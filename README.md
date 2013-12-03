gr_code
=======

code to analyse evolution of the universe

sFactor  - calculates evolution of the universe based on the 'standard  model', uses 'helper' as the helper function

sFactorCurve - solves the general Friedman equation for a universe with arbitrary curvature.  Uses 'helperCurve' as the
            the helper function
          
sFactorH - solves the Friedman equation when the equation of state parameter, w, varies as a step function.  
            Uses 'heavisideW' as the helper function.
            
sFactorA - solves the Friedman equation when the equation of state parameter, w, varies in an analytic fashion. 
            Uses 'analyticW' as the helper function.  
