function liqFraction = LiquidFraction(H,cp_solid, Tm, q)
%this function calculates the liquid function at a given temperature
%which will be used to calculate the density of the PCM

if H>=(cp_solid*Tm)
    liqFraction = 0;
elseif H>=(cp_solid*Tm +q)
    liqFraction = 1;
else 
    liqFraction = (H-cp_solid)/q;
end
end