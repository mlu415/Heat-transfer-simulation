function cp = PCMcp(T,Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition)
% calculate the cp of a PCM using an approximation method where the cp of
% the PCM ramps during a phase change between liquid and solid.


% If the PCM is fully solid
if T <= (Transition_temp - Transition_range/2)
    cp = cp_solid;
    
% If the PCM is fully liquid
elseif T >= (Transition_temp + Transition_range/2)
    cp = cp_liquid;
% If the PCM is in mushy region
else
    % Ramp functions using line equations
    if T <= Transition_temp
        cp = cp_solid + (T-(Transition_temp-Transition_range/2))*(cp_transition - cp_solid)/(Transition_range/2);
    else
        cp = cp_transition + (T-Transition_temp)*(cp_solid-cp_transition)/(0.5*Transition_range);
    end
end

end

