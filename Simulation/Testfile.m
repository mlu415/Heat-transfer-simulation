clear

%--------------------------------------------------------------------------------------------------
% Simulation parameters


nodeNumx = 10;                               % Number of nodes to break bar into
nodeNumy = 10;

%--------------------------------------------------------------------------------------------------
%SlabParameters
% Slab Diementions
L = 0.43;                                   % Length of slab
W = 0.28;                                   % Width of slab
dx = L/nodeNumx;                            % length of each node
dy = W/nodeNumy;                        % Width of slab node
dz = 0.01;                                  % Thickness of slab node
convArea = dx*dz;                           % Area for convection
condArea = dx*dy;                           % Area for conduction
nodeVol = dx*dy*dz;                         % Node volume

% Slab parameters
densitySlab=2700;                           % Density of slab
kSlab=0.205;                                  % Slab conductivity
cpSlab = 0.9;                               % Specific heat capacity kj/kgK
nodeMass = nodeVol*densitySlab;             % Node Mass
alpha=kSlab/(densitySlab*cpSlab);           % Thermal Diffusivity
initialSlabTemp=2;                          % Initial Slab temperature

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
initialPCMTemp = -4;
kPCMliq = 1.6;
kPCMsolid=0.6;

% Intialise PCM Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform

cpCurrent = PCMcp(Tp,Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
H = cpCurrent*Tp; %Calculating enthalpy at specific node and time

liqFraction = LiquidFraction(H,cp_solid,Transition_temp, cp_transition);
densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
alpha=kPCM/(densityPCM*cpCurrent);
%-----------------------------------------------------------------------------------------------------
%Recording and Simulation

simTime=400*60;  % Max simulation time
dt = min(0.5*dx^2/alpha, 0.5*dy^2/alpha);
timeSteps = round(simTime/dt);                            % Number of Time steps
fps = 0.1;
timeStepSkips = round((timeSteps/simTime)/fps);

% Intialise Slab Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialSlabTemp;                      % Setting initial slab temperature assuming uniform
%--------------------------------------------------------------------------------------------------
%AirParameters
InletAtemp=20;                              % Constant inlet air temperature degrees Celcius
h=17.2;                                     % Convective heat transfer coefficient
Ta= zeros(nodeNumx,2);               % Making empty matrix for air temperature to store values.
Ta(1,:)= InletAtemp;                        % Setting air temperature at x=
Ta(:,1)= initialSlabTemp;                        % Setting air temperature at x=0

cpAir=1.005;                                % Specific Heat capacity of air
pAir=1.205;                                 % Density of air
Velocity=2;                                 % Velocity of air
CSarea=0.0025;                              % Cross Sectional area to flow
VolumetricFlowrate=Velocity*CSarea;         % Volumetric flow rate of air
massAir=VolumetricFlowrate*pAir;                 % Mass Flow rate of air

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
q =	334; % KJ/Kg
densitySolid = 900; % Kg/m3
densityLiquid = 1000; % Kg/m3
initialPCMTemp = -4;
kPCM = 1.6;

cpCurrent = PCMcp(T,Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
H = cpCurrent*Tp; %Calculating enthalpy at specific node and time

liqFraction = LiquidFraction(H,cp_solid,Transition_temp, q);
densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;

% Intialise PCM Array
Tp= zeros(nodeNumx,nodeNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform

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
    airItr = airItr +1;
    for y = 1:nodeNumy
        for x = 1:nodeNumx                         %Change in distance
            %        %Convection
            %        if (airItr >= airSkips)
            %        Ta(x,t+1)=Ta(x,t)+(Ts(x,t)-Ta(x,t))*h*convArea*dx/(massAir*cpAir); % Air temperature profile
            %        else
            
            % Setting up next time step
            Tp(:,:,1)= Tp(:,:,2);
            Ta(:,1) = Ta(:,2);
            cpCurrent = PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
            H = cpCurrent*Tp(x,y,1); %Calculating enthalpy at specific node and time
            liqFraction = LiquidFraction(H,cp_solid,Transition_temp, cp_transition);
            densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
            kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
            
            
            alpha=kPCM/(densityPCM*cpCurrent);
            
            Stabilitycheckdx= alpha*dt/dx^2;
            Stabilitycheckdy= alpha*dt/dy^2;
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in x direction
            % Boundry conditions for slab edges
            if(x == 1)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2) + (alpha)*((Tp(x+1,y,1)-Tp(x,y,1))/dx)*dt;
            elseif(x == nodeNumx)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (alpha)*((-Tp(x,y,1)+Tp(x-1,y,1))/dx)*dt;
                Tp(x,y,2) = Tp(x,y,2)+ (alpha)*((-Tp(x,y,1)+Tp(x-1,y,1))/dx)*dt;
                % Conditions for non edge nodes
            else
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (alpha)*((Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/(dx^2))*dt;
            end
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in y direction
            % Boundary conditions for slab edge adjacent to air
            if(y==1)
                % Air temperature profile
                Ta(x+1,2)=Ta(x,1)+(Ts(x,y,1)-Ta(x,1))*h*convArea*2/(massAir*cpAir*1000); 
                % Slab Convection
                Ts(x,y,2)= Ts(x,y,2)-(massAir*cpAir*dt*(Ta(x+1,2)-Ta(x,1)))/(nodeMass*cpSlab);
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2) + (alpha)*((Ts(x,y+1,1)-Ts(x,y,1))/dy)*dt;
            % Boundry conditions for slab edge adjacent to PCM
            elseif(y == nodeNumy)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2) + (alpha)*((Tp(x,y-1,1)-Tp(x,y,1))/dy)*dt;
                % Conditions for non edge nodes
            else
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (alpha)*((Tp(x,y+1,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/(dy^2))*dt;
            end
        end
    end
    
    % Saving data into a video
    if(mod(t,timeStepSkips)==0)
%         figure(1);
%         contourf(Tp(:,:,1));
%         colorbar()
%         caxis([2, 20]);
%         drawnow;
%         F = getframe(gcf);
%         writeVideo(aviobj,F);
%         disp(100*t/timeSteps);
   end
    
    timeplot=t*dt;
    plot (timeplot,Tp(nodeNumx,nodeNumy,1),'r*-')

    
end

hold off
close(aviobj);

disp("Done");



