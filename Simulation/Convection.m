function [Q] = Convection(h, A, Tf, Ti)
%This function returns the convection with the following inputs, h 
%(Convection Coefficient), A (Area), Tf (Final Temperature) and T (Surface
%temperature)

Q = h*A*(Tf - Ti);

end

