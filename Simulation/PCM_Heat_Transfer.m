clear

%--------------------------------------------------------------------------------------------------
% Simulation parameters
simTime=2;                                 % Max simulation time
timeSteps = 100000000;                            % Number of Time steps
dt = simTime/timeSteps;                     % Time step

nodeNumx = 10;                               % Number of nodes to break bar into
nodeNumy = 10;
time = linspace(0,simTime,timeSteps);       % Array of time steps

%--------------------------------------------------------------------------------------------------
%SlabParameters
% Slab Diementions
L = 0.01;                                   % Length of slab
Thicc = 0.0001;
dx = L/nodeNumx;                            % length of each node
dy = Thicc/nodeNumy;                                % Thickness of slab
dz = 0.28;                                  % Width of slab
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


% Intialise Slab Array
Ts= zeros(nodeNumx,nodeNumy,timeSteps);               % Making empty matrix for slab to store values
Ts(:,:,1)= initialSlabTemp;                      % Setting initial slab temperature assuming uniform

%--------------------------------------------------------------------------------------------------
%AirParameters
InletAtemp=20;                              % Constant inlet air temperature degrees Celcius
h=17.2;                                     % Convective heat transfer coefficient
Ta= zeros(nodeNumx,timeSteps);               % Making empty matrix for air temperature to store values.
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

dtmin= min(0.5*dx^2/alpha, 0.5*dy^2/alpha)

airSkips = (dx/dt)/Velocity;
airItr = 0;

for t= 1:(timeSteps-1)                   %Change in Time
    airItr = airItr +1;
    for y = 1:nodeNumy
        for x = 1:nodeNumx                         %Change in distance
            %        %Convection
            %        if (airItr >= airSkips)
            %        Ta(x,t+1)=Ta(x,t)+(Ts(x,t)-Ta(x,t))*h*convArea*dx/(massAir*cpAir); % Air temperature profile
            %        else
            
%             Setting up next time step
            Ts(x,y,t+1)= Ts(x,y,t);
            
            if(x == 1)
                Ts(x,y,t+1) = Ts(x,y,t+1) + (alpha)*((Ts(x+1,y,t)-Ts(x,y,t))/dx)*dt;
            elseif(x == nodeNumx)
                Ts(x,y,t+1) = Ts(x,y,t+1)+ (alpha)*((-Ts(x,y,t)+Ts(x-1,y,t))/dx)*dt;
            else
                Ts(x,y,t+1) = Ts(x,y,t+1)+ (alpha)*((Ts(x+1,y,t)-2*Ts(x,y,t)+Ts(x-1,y,t))/(dx^2))*dt;
            end
            if(y==1)
                Ta(x+1,t+1)=Ta(x,t)+(Ts(x,y,t)-Ta(x,t))*h*convArea*2/(massAir*cpAir*1000); % Air temperature profile
                %Convection Slab
                Ts(x,y,t+1)= Ts(x,y,t+1)-(massAir*cpAir*dt*(Ta(x+1,t+1)-Ta(x,t)))/(nodeMass*cpSlab);
                Ts(x,y,t+1) = Ts(x,y,t+1) + (alpha)*((Ts(x,y+1,t)-Ts(x,y,t))/dy)*dt;
            elseif(y == nodeNumy)
                Ts(x,y,t+1) = Ts(x,y,t+1) + (alpha)*((Ts(x,y-1,t)-Ts(x,y,t))/dy)*dt;
            else
                Ts(x,y,t+1) = Ts(x,y,t+1)+ (alpha)*((Ts(x,y+1,t)-2*Ts(x,y,t)+Ts(x,y-1,t))/(dy^2))*dt;
            end
        end
    end
    
    %     if (airItr >= airSkips)
    %         airItr = 0;
    %     end
end
% TsF(:,1,:) = Ts(:,:);
% TsF(:,2,:) = Ts(:,:);
 colormap parula
% 
% aviobj = VideoWriter('example.avi');
% open(aviobj);
for t = 1:(timeSteps):50
    contourf(Ts(:,:,t));
    colorbar()
    drawnow;
%     F = getframe(gcf);
%     writeVideo(aviobj,F);
end

% aviobj = close(aviobj);

disp("Done");



