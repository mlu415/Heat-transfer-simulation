clear

%--------------------------------------------------------------------------------------------------
% Simulation parameters


nodeNumx = 5;                               % Number of nodes to break bar into
nodeNumy = 5;

%--------------------------------------------------------------------------------------------------
%SlabParameters
% Slab Diementions
% Slab Diementions
L = 0.43;                                   % Length of slab
Thicc = 0.0075;
dx = L/nodeNumx;                            % length of each node
dy = Thicc/nodeNumy;                        % Thickness of slab node
dz = 0.28;                                  % Width of slab node
convArea = dx*dz;                           % Area for convection
condArea = dx*dy;                           % Area for conduction
nodeVol = dx*dy*dz;                         % Node volume

% % Slab parameters
% densitySlab=2700;                           % Density of slab
% kSlab=0.205;                                  % Slab conductivity
% cpSlab = 0.9;                               % Specific heat capacity kj/kgK
% nodeMass = nodeVol*densitySlab;             % Node Mass
% alpha=kSlab/(densitySlab*cpSlab);           % Thermal Diffusivity
% initialSlabTemp=2;                          % Initial Slab temperature

%--------------------------------------------------------------------------------------------------
%PCM Parameters
% water parameters
Transition_temp = 0;
Transition_range = 2;

cp_liquid = 4.18; %KJ/KgK
cp_solid = 2.04; %KJ/KgK
cp_transition =	334; % KJ/Kg
densitySolid = 900; % Kg/m3
densityLiquid = 1000; % Kg/m3
initialPCMTemp = 9.999;
kPCMliq = 1.6;
kPCMsolid=0.6;

% Intialise PCM Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform
% H= zeros(nodeNumx,nodeNumy,1);
% cpCurrent = PCMcp(Tp,Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
% H(:,:,1) = cpCurrent*Tp; %Calculating enthalpy at specific node and time
liqFraction= LiquidFraction1(-1,1,Tp(:,:,:));
%Not working for some reason - liqFraction = LiquidFraction(H,cp_solid,Transition_temp, cp_transition);

densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
alpha=kPCM(1)/(densityPCM(1)*cp_solid*1000);
%-----------------------------------------------------------------------------------------------------
%Recording and Simulation

simTime=10*60;  % Max simulation time
%dt = min(0.5*dx^2/alpha, 0.5*dy^2/alpha);
dt=0.1;
timeSteps = round(simTime/dt);                            % Number of Time steps
fps = 0.1;
timeStepSkips = round((timeSteps/simTime)/fps);

% Intialise Slab Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform
%--------------------------------------------------------------------------------------------------
%AirParameters
InletAtemp=10;                              % Constant inlet air temperature degrees Celcius
h=100;                                     % Convective heat transfer coefficient
Ta= zeros(nodeNumx,2);               % Making empty matrix for air temperature to store values.
Ta(1,:)= InletAtemp;                        % Setting air temperature at x=
%Ta(:,1)= initialPCMTemp;                        % Setting air temperature at x=0

cpHTF=4.187;                                % Specific Heat capacity of air
pHTF=1000;                                 % Density of air
% Velocity of air
CSarea=pi*(0.015)^2;                              % Cross Sectional area to flow
VolumetricFlowrate=0.00005;         % Volumetric flow rate of air
Velocity=VolumetricFlowrate/CSarea;
massAir=VolumetricFlowrate*pHTF;                 % Mass Flow rate of air

%--------------------------------------------------------------------------------------------------
%row is distance where row 1 is at x=0, row 2 x=dx
%column is time where column 1 is time=0, column 2 is time =dt
%        Ta(x,t)=Ta(x-1,t)+(Ts(x,t)-Ta(x-1,t))*h*area*dx/(Ma*cpair); % Air temperature profile

% Conduction

count = 0;

dtmin= min((0.5*dx^2)/alpha, (0.5*dy^2)/alpha);

%----------------------------------------------------------------------------------------------
%PCM Parameters
% water parameters
Transition_temp = 0;
Transition_range = 2;

cp_liquid = 4.18; %KJ/KgK
cp_solid = 2.04; %KJ/KgK
cp_transition =	334; % KJ/Kg
densitySolid = 900; % Kg/m3
densityLiquid = 1000; % Kg/m3
initialPCMTemp = 9.999;
kPCMliq = 1.6;
kPCMsolid=0.6;

% Intialise PCM Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform
H = zeros (nodeNumx,nodeNumy,2);
cpCurrent = cp_solid;%PCMcp(Tp(1,1,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
H(:,:,:)= cpCurrent*1000*Tp(:,:,:);
% cpCurrent = PCMcp(Tp,Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
% H = cpCurrent*Tp; %Calculating enthalpy at specific node and time
%
% liqFraction = LiquidFraction(H,cp_solid,Transition_temp, cp_transition);
% densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
% kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;


timeplot=0;

airSkips = (dx/dt)/Velocity;
airItr = 0;
aviobj = VideoWriter('example.avi');
aviobj.FrameRate = fps;
open(aviobj);
colormap parula
figure (2)
hold on
for t= 1:(timeSteps-1)                   %Change in Time
    %airItr = airItr +1;
    %Tp(:,:,1)= Tp(:,:,2);
    %Ta(:,1)=Ta(:,2);
    
    %Ta(:,1) = Ta(:,2);
    H(:,:,1)=H(:,:,2);
    for y = 1:nodeNumy
        for x = 1:nodeNumx                         %Change in distance
            
            Tp(x,y,1)=Tp(x,y,2);
            Ta(:,1) = Ta(:,2);
            H(x,y,1)=H(x,y,2);
            
            % Setting up next time step
            
            liqFraction = LiquidFraction1(-1,1, Tp(x,y,1));
            densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
            kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
            nodeMass = nodeVol*densityPCM;
            alpha= kPCM/((H(x,y,1)/Tp(x,y,1))*densityPCM);
                         Stabilitycheckdx= alpha*dt/dx^2;
                         Stabilitycheckdy= alpha*dt/dy^2;
            
            if (Stabilitycheckdx>0.5) OR (Stabilitycheckdy>0.5)
            return
            end
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in x direction
            % Boundry conditions for slab edges
            
            if (Tp(x,y,1)-Ta(x,1))>0
                break
            end
            
            if(x == 1) && (y==1)
                
                Ta(x+1,2)=Ta(x,2)+(Tp(x,y,1)-Ta(x,1))*h*convArea*2/(massAir*cpHTF*1000);
                
                % PCM Convection
                H(x,y,1)=H(x,y,1)-(massAir*cpHTF*dt*(Ta(x+1,2)-Ta(x,2))/nodeMass);
                
                Tp(x,y,1)=PCMcp2(H(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+2,y,1)-2*Tp(x+1,y,1)+Tp(x,y,1))/dx^2)+((Tp(x,y+2,1)-2*Tp(x,y+1,1)+Tp(x,y,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
            elseif(x == nodeNumx) && (y==nodeNumy)
                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x-2,y,1)-2*Tp(x-1,y,1)+Tp(x,y,1))/dx^2)+((Tp(x,y-2,1)-2*Tp(x,y-1,1)+Tp(x,y,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
                
            elseif(x ~= 1) && (x~=nodeNumx) && (y==nodeNumy)
                                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/dx^2)+((Tp(x,y-2,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
                
            elseif(x == nodeNumx) && (y == 1)
                
                Ta(x+1,2)=Ta(x,2)+(Tp(x,y,1)-Ta(x,2))*h*convArea*2/(massAir*cpHTF*1000);
                
                % PCM Convection
                H(x,y,1)=H(x,y,1)-(massAir*cpHTF*dt*(Ta(x+1,2)-Ta(x,2))/nodeMass);
                
                Tp(x,y,1)=PCMcp2(H(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x-2,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/dx^2)+((Tp(x,y+2,1)-2*Tp(x,y+1,1)+Tp(x,y,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                
                
            elseif(x == 1) && (y==nodeNumy)
               
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+2,y,1)-2*Tp(x+1,y,1)+Tp(x,y,1))/dx^2)+((Tp(x,y-2,1)-2*Tp(x,y-1,1)+Tp(x,y,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
                
            elseif(y==1) && (x~=1) && (x~=nodeNumx)
                
                Ta(x+1,2)=Ta(x,2)+(Tp(x,y,1)-Ta(x,2))*h*convArea*2/(massAir*cpHTF*1000);
                
                % PCM Convection
                H(x,y,1)=H(x,y,1)-(massAir*cpHTF*dt*(Ta(x+1,2)-Ta(x,2))/nodeMass);
                
                Tp(x,y,1)=PCMcp2(H(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/dx^2)+((Tp(x,y+2,1)-2*Tp(x,y+1,1)+Tp(x,y,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                
            elseif(x==1) && (y~=1) && (y~=nodeNumy)
                                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+2,y,1)-2*Tp(x+1,y,1)+Tp(x,y,1))/dx^2)+((Tp(x,y+1,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
                
            elseif(x==nodeNumx) && (y~=1) && (y~=nodeNumy)
                
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x-2,y,1)-2*Tp(x-1,y,1)+Tp(x,y,1))/dx^2)+((Tp(x,y+1,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
                
            else
               
                % PCM Conduction
                H(x,y,2)=(kPCM*dt/densityPCM)*(((Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/dx^2)+((Tp(x,y+1,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/dy^2))+ H(x,y,1);
                
                Tp(x,y,2)=PCMcp2(H(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid);
                
                Tp(x,y,1)=Tp(x,y,2);
            end
            
        end
        
    end
    
    %Saving data into a video
%     if(mod(t,timeStepSkips)==0)
%         figure(1);
%         contourf(Tp(:,:,1));
%         colorbar()
%         caxis([1, 20]);
%         drawnow;
%         F = getframe(gcf);
%         writeVideo(aviobj,F);
%         disp(100*t/timeSteps);
%     end
    
        timeplot=t*dt;
        plot (timeplot,Tp(1,1,1),'r*-')
    
    
end

hold off
close(aviobj);

disp("Done");



