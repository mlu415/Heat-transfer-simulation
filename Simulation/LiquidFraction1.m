function liqFraction1 = LiquidFraction1(Tps,Tpl,T)
%this function calculates the liquid function at a given temperature
%which will be used to calculate the density of the PCM

if T<Tps
    liqFraction1 = 0;
elseif T>=Tpl
    liqFraction1 = 1;
else 
    liqFraction1 = (T-Tps)/(Tpl-Tps);
end
end