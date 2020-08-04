function T = PCMcp2(H,Transition_temp,Transition_range,cp_liquid,cp_solid)%,cp_transition)
%function T = PCMcp2(H,Transition_temp,Transition_range,cp_liquid,cp_solid)%,cp_transition)
% calculate the cp of a PCM using an approximation method where the cp of
% the PCM ramps during a phase change between liquid and solid.


% If the PCM is fully solid
if H/(1000*cp_solid) < 0%(Transition_temp - Transition_range/2)
    T=H/(1000*cp_solid);
    
    
    % If the PCM is fully liquid
elseif H/(1000*cp_liquid) > 0%(Transition_temp + Transition_range/2)
    T=H/(1000*cp_liquid);
    
    % If the PCM is in mushy region
else
    % Ramp functions using line equations (interpolation)
    T=0 + 0*(Transition_temp+Transition_range);
    %cp = cp_solid + (T-(Transition_temp-Transition_range/2))*(cp_transition - cp_solid)/(Transition_range/2);
end


end

