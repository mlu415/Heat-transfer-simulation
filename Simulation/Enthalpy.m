clear

%--------------------------------------------------------------------------------------------------
% Simulation parameters

slabNumx = 10;                               % Number of nodes to break bar into
slabNumy = 10;

%--------------------------------------------------------------------------------------------------
%HTFParameters
InletHtfTemp=10 + 273;                              % Constant inlet air temperature degrees Celcius
h=17.2;                                     % Convective heat transfer coefficient
 Ta= zeros(slabNumx,1);               % Making empty matrix for air temperature to store values.
 Ta(:)= InletHtfTemp;                        % Setting air temperature at x=

cpHTF=1.005;                                % Specific Heat capacity of air
pHTF=1.205;                                 % Density of air
Velocity=2;                                 % Velocity of air
viscosityHTF = 1.511e-5;                    % Kinematic viscosity
kHTF = 0.0257;                              % Thermal conductivity of air(W/m.K)
CSarea=0.0025;                              % Cross Sectional area to flow
VolumetricFlowrate=Velocity*CSarea;         % Volumetric flow rate of air
massHTF=VolumetricFlowrate*pHTF;            % Mass Flow rate of air

%--------------------------------------------------------------------------------------------------
%SlabParameters
% Slab Dimensions
L = 0.43;                                   % Length of slab
thickSlab = 0.01;
dx = L/slabNumx;                            % length of each node
dy = thickSlab/slabNumy;                        % Thickness of slab node
dz = 0.28;                                  % Width of slab node
nodeVol = dx*dy*dz;                         % Node volume

% Slab parameters
densitySlab=2700;                           % Density of slab
kSlab=0.205;                                  % Slab conductivity
cpSlab = 0.9;                               % Specific heat capacity kj/kgK
nodeMass = nodeVol*densitySlab;             % Node Mass
alpha=kSlab/(densitySlab*cpSlab);           % Thermal Diffusivity
initialSlabTemp=10 + 273;                          % Initial Slab temperature
gap = 0.005;                                %gap between slabs
nSlab = 30;
velocityAir = (Velocity*pi*(0.1/2)^2)/L*gap*nSlab;
charLength=(2*dz*gap)/(dz+gap);

Re=(velocityAir*charLength)/viscosityHTF;
Pr=(cpHTF*viscosityHTF*pHTF*1000)/kHTF;      % Prandtl number of the air
Nu=Nusselts(dz, gap,L, charLength,Re,Pr);
ha = (Nu*kHTF/charLength);

if (Re <= 2300)                               % Darcy friction factor for parallel plates
    f=64/(2/3*Re);
elseif (Re > 2300) && (Re < 1e6)
    f=(0.76*log(Re)-1.64)^(-2);
end



%-----------------------------------------------------------------------------------------------------
%Recording and Simulation

simTime=10*60;                                 % Max simulation time
dt = min(0.5*dx^2/alpha, 0.5*dy^2/alpha)/10;
timeSteps = round(simTime/dt);                            % Number of Time steps
fps = 100;
timeStepSkips = round((timeSteps/simTime)/fps);

% Intialise Slab Array
Ts= zeros(slabNumx,slabNumy,2);               % Making empty matrix for slab to store values
Ts(:,:,:)= initialSlabTemp;                      % Setting initial slab temperature assuming uniform
%Ta(:,1)= initialSlabTemp;                        % Setting air temperature at x=0

%----------------------------------------------------------------------------------------------
%PCM Parameters
AvgPCMtemp = zeros(timeSteps,1);
PointPCMtemp = zeros(timeSteps,1);
% water parameters
Transition_temp = 0 + 273;
Transition_range = 2;

cp_liquid = 4.18; %KJ/KgK
cp_solid = 2.04; %KJ/KgK
cp_transition =	334; % KJ/Kg

latentH = 334;

kPCMliq = 1.6;
kPCMsolid=0.6;
hPCM = 50;
densitySolid = 900; % Kg/m3
densityLiquid = 1000; % Kg/m3
initialPCMTemp = 10+ 273;
kPCM = 1.6;
PCMthickness = 0.01; % 1cm
PCMNumy = 10;
PCMdy = PCMthickness/PCMNumy;

H_PCM_solid = cp_solid*(Transition_temp-Transition_range/2);                           % Enthalpy of PCM in solid phase (kJ/kg)
H_PCM_liquid = (cp_solid*(Transition_temp))+cp_transition+Transition_range/2*cp_liquid;                % Enthalpy of PCM in liquid phase (kJ/kg)




% Intialise PCM Array
Tp= zeros(slabNumx,PCMNumy,2);               % Making empty matrix for slab to store values
Tp(:,:,:)= initialPCMTemp;                      % Setting initial slab temperature assuming uniform
Tsol = -1;
Tliq = 1;

Hp = zeros(slabNumx,PCMNumy,2);
H(:,:,:)= cp_solid*Tp(:,:,:);

% if (Tp(:,:) < Transition_temp-(Transition_range/2))
%     Hp(:,:)=cp_solid.*Tp(:,:);
% elseif (Transition_temp-(Transition_range/2) <= Tp(x,y)) && ( Tp(x,y)<= Transition_temp+(Transition_range/2))
%     Hp(:,:)=cp_solid*(Transition_temp-Transition_range/2)+((cp_solid+cp_liquid)/2+(cp_transition/(Transition_range))).*(Tp(:,:)-Transition_temp+Transition_range/2);
% elseif (Tp(:,:) > Transition_temp+Transition_range/2)
%     Hp(:,:)=cp_liquid.*Tp(:,:)+(cp_solid-cp_liquid)*Transition_temp+cp_transition;
% end



airSkips = (dx/dt)/Velocity;
airItr = 0;
aviobj = VideoWriter('example.avi');
aviobj.FrameRate = fps;
open(aviobj);
colormap parula




for t= 1:(timeSteps-1)                   %Change in Time
    Ts(:,:,1)= Ts(:,:,2);
    Ta(:) = Ta(:);
    Tp(:,:,1)= Tp(:,:,2);
    H(:,:,1)=H(:,:,2);
    airItr = airItr +1;
    
    
    for x = 1:slabNumx                         %Change in distance
        for y = 1:PCMNumy
            
            liqFraction = LiquidFractionE(Tsol,Tliq,Tp(x,y,1));
            densityPCM = liqFraction*densityLiquid + (1-liqFraction)*densitySolid;
            kPCM = liqFraction*kPCMliq + (1-liqFraction)*kPCMsolid;
            pcmAlpha = kPCM/((H(x,y,1)/Tp(x,y,1))*densityPCM);
            massPCM = densityPCM*dx*PCMdy*dz;
            cpPCM = cp_liquid*liqFraction+(1-liqFraction)*cp_solid;
            
            Stabilitycheckdx= pcmAlpha*dt/dx^2;
            Stabilitycheckdy= pcmAlpha*dt/PCMdy^2;
            if(Stabilitycheckdx >= 0.5 || Stabilitycheckdy >= 0.5)
                fprintf('Error')
            end
            
            % -------------------------------------------------------------
            % Calulating temperature change from heat flow
            % Boundry conditions for PCM edges
            if (y == 1)
                %                 Tp(x,y)=Tp(x,y)+3.6*dt/((densitySlab*cpSlab*thickSlab)+(densityPCM(x,y).*cpPCM(x,y)*dz))*(0.85*ha*(Ta(x)-Tp(x,y))-kPCM(x,y)/dz*(Tp(x,y)-Tp(x,y+1)));
                
                if (x==1)
                    Tp(x,y,2)=Tp(x,y,1)+3.6*dt/(densitySlab*cpSlab*thickSlab+(densityPCM*cpPCM*dz))*(0.85*ha*(Ta(x)-Tp(x,y,1))-kPCM/dz*(Tp(x,y,1)-Tp(x,y+1,1))+kSlab*thickSlab/dx^2*(Tp(x+1,y,1)-Tp(x,y,1)));
                elseif (2 <= x) && (x < slabNumx)
                    Tp(x,y,2)=Tp(x,y,1)+3.6*dt/(densitySlab*cpSlab*thickSlab   +    (densityPCM*cpPCM*dz)   )*(  0.85*ha*(Ta(x)-Tp(x,y,1))  -  kPCM/dz*(Tp(x,y,1)-Tp(x,y+1,1))  +   kSlab*thickSlab/dx^2*(Tp(x+1,y,1)-2*Tp(x,y,1)+Tp(x-1,y,1)));
                elseif (x==slabNumx)
                    Tp(x,y,2)=Tp(x,y,1)+3.6*dt/(densitySlab*cpSlab*thickSlab+(densityPCM*cpPCM*dz))*(0.85*ha*(Ta(x)-Tp(x,y,1))-kPCM/dz*(Tp(x,y,1)-Tp(x,y+1,1))+kSlab*thickSlab/dx^2*(-Tp(x,y,1)+Tp(x-1,y,1)));
                end
                
                
            elseif ((2 <= y) && (y < PCMNumy))
                
                if (x==1)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y+1,1)-2.*Tp(x,y,1)+Tp(x,y-1,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(Tp(x+1,y,1)-Tp(x,y,1));
                elseif (2 <= x) && (x < slabNumx)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y+1,1)-2.*Tp(x,y,1)+Tp(x,y-1,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(Tp(x+1,y,1)-2.*Tp(x,y,1)+Tp(x-1,y,1));
                elseif (x==slabNumx)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y+1,1)-2.*Tp(x,y,1)+Tp(x,y-1,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(-Tp(x,y,1)+Tp(x-1,y,1));
                end
                
                
            elseif (y == PCMNumy)
                if (x==1)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y-1,1)-Tp(x,y,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(Tp(x+1,y,1)-Tp(x,y,1));
                elseif (2 <= x) && (x < slabNumx)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y-1,1)-Tp(x,y,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(Tp(x+1,y,1)-2.*Tp(x,y,1)+Tp(x-1,y,1));
                elseif (x==slabNumx)
                    Hp(x,y)=Hp(x,y)+((kPCM*dt*3.6)./(densityPCM*(dz^2))).*(Tp(x,y-1,1)-Tp(x,y,1))+((kPCM*dt*3.6)./(densityPCM*(dx^2))).*(-Tp(x,y,1)+Tp(x-1,y,1));
                end       
            end
            
            if (Hp(x,y) < H_PCM_solid)
                Tp(x,y,1)=Hp(x,y)/cp_solid;
            elseif (H_PCM_solid <= Hp(x,y)) && (Hp(x,y) < H_PCM_liquid)
                Tp(x,y,1)=(Hp(x,y)-cp_solid*(Transition_temp-Transition_range/2))/((cp_solid+cp_liquid)/2+latentH/(Transition_range))+Transition_temp-(Transition_range/2);
            elseif (Hp(x,y)>= H_PCM_liquid)
                Tp(x,y,1)=(Hp(x,y)-(cp_solid)*Transition_temp-latentH)/cp_liquid+Transition_temp;
            end
            
        end
    end
    
    % Saving data into a video
    if(mod(t,timeStepSkips)==0)
%                     hold on
%                     contourf(flipud(Tp(:,:,1).'));
%                     colorbar()
%                     caxis([-100, 20]);
%                     drawnow;
%                     F = getframe(gcf);
%                     writeVideo(aviobj,F);
%         disp(100*t/timeSteps);
%         hold off
    end
    AvgPCMtemp(t) = mean(Tp(:,:,1), "all");
    PointPCMtemp(t) = Tp(10,2,1);
end

figure
hold on
plot((1:(timeSteps-1))*dt/60, AvgPCMtemp(1:end-1));
plot((1:(timeSteps-1))*dt/60, PointPCMtemp(1:end-1));
xlabel('Time (Mins)');
ylabel('Temperature (C)')
hold off
close(aviobj);

disp("Done");



