function Nu=Nusselts(dz, gap,length, dH, Re, Pr)

if (Re < 2300) && (dz/gap > 8)
    Nu=7.541+(0.03*(dH/length)*Re*Pr)/(1+0.016*((dH/length)*Re*Pr).^(2/3));        % Nusselt number
elseif (Re >= 2300) && (Re < 10e6)
    f=(0.76*log(Re)-1.64)^(-2);
    Nu=((f/8)*(Re-1000)*Pr)/(1+12.7*(f/8)^0.5*(Pr^(2/3)-1));
end

