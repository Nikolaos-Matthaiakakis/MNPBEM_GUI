% *********************************************************
% ITO optical model 
% *********************************************************
%Matthaiakakis, Nikolaos      (2017)     
%Dynamic modulation of plasmon excitations in monolayer graphene.  
%University of Southampton, Doctoral Thesis, 231pp. 
function ITO_m( )
% *********************************************************
% Input values
% *********************************************************
%T=300; % Temperature K
epsilon_inf=3.9; %infinate frequency permittivity
omega_p=0.5047; %wavelength range in eV
gama_=20; %wavelength range in um
ene2=0.6:0.005:1.2; %wavelength range in um
%la=0.6:0.005:1.2; %wavelength range in um
%------ 
%n_cr=0.5;%:1:5; %carrier density *10^20 cm^3
%------ Voltage input
%V=0;%:0.005:5; % Voltage V
%d_oxide=10*10^-9; %m thickness of diel
%epsilon_oxide=25; %permittivity of dielectric
%epsilon_ITO=3.34; %permittivity of dITO
%A_r=1; %contact area, leave at 1 for capacitance in F/cm^2. otherwise use value in cm^2

set_=3;%calculates from voltage, 2 from carrier conc, other from omega_p
% *********************************************************
% Constants
% *********************************************************
h_bar=1.05*10^-34; %J.s
h_plank=6.626*10^-34; %J.s
h_plankeV=4.135667662*10^-15; %eV
epsilon_0=8.85*10^-12; %F/m permittivity of vacuum
e=1.6*10^-19; %C electron charge
c = 2.99792458*10^8; %m/s speed of light
me_s=9.109*10^-31; %electron mass in free space kg
me=0.45*me_s; %electron mass in ITO kg
% *********************************************************
% Carrier concentration in ITO
% *********************************************************
n_ct=n_cr*10^20*(2*pi)^2;%1/cm^3
n_c=(n_ct/10^-6);%convert to 1/m^3
n_c2=(n_c/(2*pi)^2+(epsilon_0*epsilon_oxide*V/(e*d_oxide^2)))*(2*pi)^2;
% *********************************************************
% Layer thickness
% *********************************************************
tTF=((pi^4./(3*n_c)).^(1/6)).*(epsilon_0*epsilon_ITO*h_plank^2/(4*me*e^2*pi^2))^(1/2);
% *********************************************************
% Permittivity and Refractive index
% *********************************************************
if set_==1
    omegap=sqrt(n_c2*e^2/((me*epsilon_0)));
elseif set_==2
    omegap=sqrt(n_c*e^2/((me*epsilon_0)));
else
    lambda_op=omega_p*10^-6; %m
    omegap=2*pi*c./lambda_op;
end
lambda=la*10^-6; %m
omega_=2*pi*c./lambda;
lambda_G=gama_*10^-6; %m
gama=2*pi*c./lambda_G;
for i=1:1:length(V)
    perm(i,:)=epsilon_inf-omegap(i).^2./(omega_.^2+1i.*gama.*omega_);
end
% *********************************************************
%Export data to text
% *********************************************************
n=sqrt(perm);
nr=real(n);
ni=imag(n);
ITO_nk=table(ene2',nr',ni');
writetable(ITO_nk,'MNPBEM17/Material/@epstable/ITO.dat','Delimiter',' ','WriteVariableNames',0)
