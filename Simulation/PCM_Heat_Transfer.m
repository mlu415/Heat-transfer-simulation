clear

%--------------------------------------------------------------------------------------------------
% Simulation parameters

slabNumx = 10;                               % Number of nodes to break bar into
slabNumy = 10;

%--------------------------------------------------------------------------------------------------
%SlabParameters
% Slab Dimensions
L = 0.43;                                   % Length of slab m
Thick = 0.015;
dx = L/slabNumx;                            % length of each node
dy = Thick/slabNumy;                        % Thickness of slab node
dz = 0.28;                                  % Width of slab node
nodeVol = dx*dy*dz;                         % Node volume

% Slab parameters
densitySlab=913;                           % Density of slab
kSlab=0.3;                                  % Slab conductivity
cpSlab = 1.67;                               % Specific heat capacity kj/kgK
nodeMass = nodeVol*densitySlab;             % Node Mass
alpha=kSlab/(densitySlab*cpSlab);           % Thermal Diffusivity
initialSlabTemp=-10;                          % Initial Slab temperature

%-----------------------------------------------------------------------------------------------------
%Recording and Simulation

simTime=30*60;                                 % Max simulation time
dt = min(0.5*dx^2/alpha, 0.5*dy^2/alpha)/10;
timeSteps = round(simTime/dt);                            % Number of Time steps
fps = 0.1;
timeStepSkips = round((timeSteps/simTime)/fps);

% Intialise Slab Array
Ts= zeros(slabNumx,slabNumy,2);               % Making empty matrix for slab to store values
Ts(:,:,:)= initialSlabTemp;                      % Setting initial slab temperature assuming uniform

%--------------------------------------------------------------------------------------------------
%HTFParameters
InletHTFtemp=5;                              % Constant inlet air temperature degrees Celcius
h=10600;                                     % Convective heat transfer coefficient 10622.6
Ta= zeros(slabNumx,2);               % Making empty matrix for air temperature to store values.
Ta(1,:)= InletHTFtemp;                        % Setting air temperature at x=
Ta(:,1)= initialSlabTemp;                        % Setting air temperature at x=0

cpHTF=3.72;                                % Specific Heat capacity of air 3.72
pHTF=1060;                                 % Density of air 1045
Velocity=0.04;                                 % Velocity of air
CSarea=0.0025;                              % Cross Sectional area to flow
VolumetricFlowrate=Velocity*CSarea;         % Volumetric flow rate of air
massAir=VolumetricFlowrate*pHTF;                 % Mass Flow rate of air

%----------------------------------------------------------------------------------------------
%PCM Parameters
AvgPCMtemp = zeros(timeSteps,1);
PointPCMtemp = zeros(timeSteps,1);
% water parameters
Transition_temp = 0;
Transition_range = 2;

cp_liquid = 4.18; %KJ/KgK
cp_solid = 2.04; %KJ/KgK
cp_transition =	334; % KJ/Kg

kPCMliq = 0.6;
kPCMsolid=0.6;
hPCM = 50;
densitySolid = 900; % Kg/m3
densityLiquid = 1000; % Kg/m3
initialPCMTemp = -10;
PCMthickness = 0.02; % 1cm
PCMNumy = 10;
PCMdy = PCMthickness/PCMNumy;

% Intialise PCM Array
Tp= zeros(slabNumx,PCMNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform

airSkips = (dx/dt)/Velocity;
airItr = 0;
aviobj = VideoWriter('example.avi');
aviobj.FrameRate = fps;
open(aviobj);
colormap parula
for t= 1:(timeSteps-1)                   %Change in Time
    Ts(:,:,1)= Ts(:,:,2);
    Ta(:,1) = Ta(:,2);
    Tp(:,:,1)= Tp(:,:,2);
    airItr = airItr +1;
    Ta(1,1) = InletHTFtemp;
    for x = 1:slabNumx                         %Change in distance
        for y = 1:PCMNumy
            
            cpPCM = PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
            H = PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition)*Tp(x,y,1); %Calculating enthalpy at specific node and time
            liqFraction = LiquidFraction(H,cp_solid,Transition_temp, cp_transition);
            densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
            kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
            pcmAlpha = kPCM/(densityPCM*cpPCM);
            massPCM = densityPCM*dx*PCMdy*dz;
            
            Stabilitycheckdx= pcmAlpha*dt/dx^2;
            Stabilitycheckdy= pcmAlpha*dt/PCMdy^2;
            if(Stabilitycheckdx >= 0.5 || Stabilitycheckdy >= 0.5)
                fprintf('Error')
            end
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in x direction
            % Boundry conditions for PCM edges
            if(x == 1)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2) + (pcmAlpha)*((Tp(x+1,y,1)-Tp(x,y,1))/dx)*dt;
            elseif(x == slabNumx)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (pcmAlpha)*((-Tp(x,y,1)+Tp(x-1,y,1))/dx)*dt;
                % Conditions for non edge nodes
            else
                % Slab Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (pcmAlpha)*((Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1))/(dx^2))*dt;
            end
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in y direction
            % Boundry conditions for PCM edge adjacent to Slab
            if(y==1)
                % PCM/Slab Convection Boundary Conditions
                Tp(x,y,2)= Tp(x,y,2)+(Ts(x,slabNumy,1)-Tp(x,y,1))*hPCM*dx*dz*2/(massPCM*cpPCM*1000);
                Ts(x,slabNumy,2)= Ts(x,y,2)-(massPCM*cpPCM*dt*(Tp(x,y,2)-Tp(x,y,1)))/(nodeMass*cpSlab);
                
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2) + (pcmAlpha)*((Tp(x,y+1,1)-Tp(x,y,1))/PCMdy)*dt;
                % Boundry conditions for slab edge adjacent to PCM
            elseif(y == PCMNumy)
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2) + (pcmAlpha)*((Tp(x,y-1,1)-Tp(x,y,1))/PCMdy)*dt;
                % Conditions for non edge nodes
            else
                % PCM Conduction
                Tp(x,y,2) = Tp(x,y,2)+ (pcmAlpha)*((Tp(x,y+1,1)-2*Tp(x,y,1)+Tp(x,y-1,1))/(PCMdy^2))*dt;
            end
        end
        
        for y = 1:slabNumy
            
            % Setting up next time step
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in x direction
            % Boundry conditions for slab edges
            if(x == 1)
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2) + (alpha)*((Ts(x+1,y,1)-Ts(x,y,1))/dx)*dt;
            elseif(x == slabNumx)
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2)+ (alpha)*((-Ts(x,y,1)+Ts(x-1,y,1))/dx)*dt;
                % Conditions for non edge nodes
            else
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2)+ (alpha)*((Ts(x+1,y,1)-2*Ts(x,y,1)+Ts(x-1,y,1))/(dx^2))*dt;
            end
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow in y direction
            % Boundry conditions for slab edge adjacent to air
            if(y==1)
                % Air temperature profile
                Ta(x+1,2)=Ta(x,1)+(Ts(x,y,1)-Ta(x,1))*h*dx*dz*2/(massAir*cpHTF*1000);
                % Slab Convection
                Ts(x,y,2)= Ts(x,y,2)-(massAir*cpHTF*dt*(Ta(x+1,2)-Ta(x,1)))/(nodeMass*cpSlab);
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2) + (alpha)*((Ts(x,y+1,1)-Ts(x,y,1))/dy)*dt;
                % Boundry conditions for slab edge adjacent to PCM
            elseif(y == slabNumy)
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2) + (alpha)*((Ts(x,y-1,1)-Ts(x,y,1))/dy)*dt;
                % Conditions for non edge nodes
            else
                % Slab Conduction
                Ts(x,y,2) = Ts(x,y,2)+ (alpha)*((Ts(x,y+1,1)-2*Ts(x,y,1)+Ts(x,y-1,1))/(dy^2))*dt;
            end
        Hnew = Tp(x,y,1)*PCMcp(Tp(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
        %testyboi = PCMcp(Tp(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition)-PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition)
        
%           Hdiff = (Tp(x,y,2)-Tp(x,y,1))*(PCMcp(Tp(x,y,2),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition)-PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition));
%           Tp(x,y,2) = Tp(x,y,2) - 1*Hdiff/PCMcp(Tp(x,y,1),Transition_temp,Transition_range,cp_liquid,cp_solid,cp_transition);
        end
        
        
    end
    
    % Saving data into a video
        if(mod(t,timeStepSkips)==0)
            contourf(flipud(cat(2,Ts(:,:,1),Tp(:,:,1)).'));
            colorbar()
            caxis([-5, 20]);
            drawnow;
            F = getframe(gcf);
            writeVideo(aviobj,F);
            disp(100*t/timeSteps);
        end
        AvgPCMtemp(t) = mean(Tp(:,:,1), "all");
        PointPCMtemp(t) = Tp(10,10,1);
end

% figure
% hold on
% plot((1:(timeSteps-1))*dt/60, AvgPCMtemp(1:end-1));
% plot((1:(timeSteps-1))*dt/60, PointPCMtemp(1:end-1));
% xlabel('Time (Mins)');
% ylabel('Temperature (C)')

close(aviobj);

disp("Done");



