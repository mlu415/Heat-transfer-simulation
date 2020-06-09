clc
clear

%System Conditions

%--------------------------------------------------------------------------------------------------
%AirParameters
InletAtemp=20;                              %Constant inlet air temperature degrees Celcius
h=17.2;                                     %Convective heat transfer coefficient
Ta= zeros(10000,10000);                     %Making empty matrix for air temperature to store values.
Ta(1,1)= InletAtemp;                        %Setting air temperature at x=0

cpair=1.005;                                %Specific Heat capacity of air
pair=1.205;                                 %Density of air
Velocity=2;                                 %Velocity of air
CSarea=0.0025;                              %Cross Sectional area to flow
VolumetricFlowrate=Velocity*CSarea;         %Volumetric flow rate of air
Ma=VolumetricFlowrate*pair;                 %Mass Flow rate of air
%--------------------------------------------------------------------------------------------------

%SlabParameters
InitialStemp=2;                             %Initial Slab temperature
Ts= zeros(10000,10000);                  %Making empty matrix for slab to store values
Ts(:,:)= InitialStemp;                   %Setting initial slab temperature assuming uniform

p=2700;                                     %Density of slab
kslab=205;                                  %Slab conductivity
cp=100;                                     %Specific heat capacity
t=0.0001;                                   %Thickness of slab
w=0.28;                                     %Width of slab
L=0.43;                                     %Length of slab in x-direction
dx=0.0001;                                  %steps in x direction
I=L/dx;                                     %Number of nodes in x-direction
dt=0.0001;                                     %Fixed Time steps
%N=L/dt;                                    %Number of nodes for time
area=w*dx;                                  %Differential area
mass=p*area*t;                              %Differential mass
tfinal=40;                                  %Simulation time
time=0:dt:tfinal;
Length=0:dx:L;
Lastnode=length(Length);
%--------------------------------------------------------------------------------------------------
%row is distance where row 1 is at x=0, row 2 x=dx
%column is time where column 1 is time=0, column 2 is time =dt

for j=1:length(time)                        %Change in Time
    Ts(1,j)=h*dx*(Ta(1,j)-Ts(1,j))/kslab+Ts(1,j); %bc1
    for i=1:Lastnode-1                              %Change in distance
        Ta(i+1,j)=Ta(i,j)+(Ts(i,j)-Ta(i,j))*h*w*dx/(Ma*cpair);
    end
    Ta(Lastnode,j)=Ta(i,j)+(Ts(i,j)-Ta(i,j))*h*w*dx/(Ma*cpair);
    Ts(Lastnode,j)=h*dx*(Ta(Lastnode,j)-Ts(Lastnode,j))/kslab+Ts(Lastnode-1,j);
   
    for i=2:Lastnode-1                              %Change in distance
       Ts(i,j+1)=(kslab*dt/(p*cp*dx^2))*(Ts(i+1,j)-2*Ts(i,j)+Ts(i-1,j));
    end
    
end

