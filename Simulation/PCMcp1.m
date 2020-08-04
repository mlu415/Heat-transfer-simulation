function T = PCMcp1(H,Transition_temp,Transition_range,cp_liquid,cp_solid)%,cp_transition)
% calculate the cp of a PCM using an approximation method where the cp of
% the PCM ramps during a phase change between liquid and solid.


% If the PCM is fully solid
if H/cp_solid <= (Transition_temp - Transition_range/2)
    T=H/cp_solid;
    
    
% If the PCM is fully liquid
elseif H/cp_liquid >= (Transition_temp + Transition_range/2)
    T=H/cp_liquid;
    
% If the PCM is in mushy region
else
    % Ramp functions using line equations (interpolation)
   syms x
    if H<= cp_solid*-1+167
        y=162.96*x+165;
        out=int(y);
        areaoftriangle=H-(cp_solid*-1);
        eqn=out==areaoftriangle;
        subs(out,0);
        solx=solve(eqn,x);
        x=-1/solx;
        T=x(1);
        %cp = cp_solid + (T-(Transition_temp-Transition_range/2))*(cp_transition - cp_solid)/(Transition_range/2);
    else
        y=-160.813*x+165;
        out=int(y);
        areaoftriangle=H-165;
        eqn=out==areaoftriangle;
        subs(out,0);
        solx=solve(eqn,x);
        x=-1/solx;
        T=x;
        %cp = cp_transition + (T-Transition_temp)*(cp_solid-cp_transition)/(0.5*Transition_range);
        
    end
end

end

