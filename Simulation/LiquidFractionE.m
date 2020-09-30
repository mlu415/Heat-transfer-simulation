function liqFraction = LiquidFractionE(Tsol,Tliq,T)
%this function calculates the liquid function at a given temperature
%which will be used to calculate the density of the PCM

if T<Tsol
    liqFraction = 0;
elseif T>=(Tliq)
    liqFraction = 1;
else 
    liqFraction = (T-Tsol)/(Tliq-Tsol);
end
end