% *********************************************************
% Graphene optical model 
% *********************************************************
%G. W. Hanson, "Dyadic Green’s Functions and Guided Surface Waves 
%for a Surface Conductivity Model of Graphene," J. Appl. Phys. 103, 
%064302; DOI:10.1063/1.2891452 (2007). 
%T. Stauber, N. M. R. Peres, and a. K. Geim, "Optical conductivity of 
%raphene in the visible region of the spectrum," Phys. Rev. B - Condens
%. Matter Mater. Phys. 78, 1–8 (2008).
function graphene_m( )%ene2, d_length2, Ef_1, mobil, T,  Drude)
% *********************************************************
% Code inputs
% *********************************************************
% Load input variables   
mat = dir('tmp/*.mat');
for q = 1:length(mat) 
    load(sprintf(['tmp/' mat(q).name])); 
end 
% *********************************************************
% Input values
% *********************************************************
mass=1.6*10^-20;
%----------------------------------------------------------
%Main program**********************************************
%----------------------------------------------------------

% *********************************************************
% Constants and variables
% *********************************************************
h_plank=6.582*10^-16; %eV.s
%epsilon_0=8.85*10^-14; %F/cm permittivity of vacuum
e=1.6*10^-19; %C electron charge
k=8.6*10^-5; %eVK^-1 Boltzmann constant
c = 2.99792458*10^8; %m/s speed of light
uf=10^8; %cm/s fermi velocity
t=2.7; %eV hopping parameter
dgr=edge_profile*10^-9; %m thickness of the graphene layer
%epsilon_air=1; % relative permittivity of air
%mv=4*pi*10^-7; %megnetic permeabillity of vacuum
%Lorentz parameters----
FL=2.234; %oscillation strength f /5
EL=4.61; %position E0
GL=0.6*2; %width gama *1.4
%mass=1.1*10^-30;%kg
ene2= linspace( minene2, maxene2, stepene2 );
if length(ene2)==1
    ene2=[ene2 ene2+0.01];
end
lambda=((4.136*10^-15)*3*10^8)./ene2; %m
omega=2*pi*c./lambda; %hz
% *********************************************************
% Fermi Level in Graphene
% *********************************************************
Ef=Ef_1; %eV
Ef1(length(Ef),length(lambda))=0;
for i=1:1:length(lambda)
    Ef1(:,i)=Ef';
end
%tmob=mobil*Ef/(e*uf^2)*(1.602177*10^-19); % plasmon lifetime s
Gama=e/(mobil*mass); %1*10^-3;% eV damping rate
if length(Ef)==1
    ngama=(Ef/(h_plank*uf))^2/pi;
    Gama=uf/(mobil*sqrt(pi*ngama));%*e eV damping rate
end
% *********************************************************
% Refractive index of graphene
% *********************************************************
cond_gr=zeros(length(Ef),length(omega)); %pre-allocation of matrices
cond_gi=cond_gr;
cond_g=cond_gr;
perm=cond_gr;
ng=cond_gr;
permlorenrz=cond_gr;
for ji=1:1:max(length(Ef)) % calculation of the optical conductivity, permitivvity and refractive index
    if Drude==0
        cond_gr(ji,:)=((e^2/(4*h_plank)*(1+(h_plank*omega).^2./(36*t^2))/2.*  (tanh((h_plank*omega+2*Ef(ji))/(4*k*T))+tanh((h_plank*omega-2*Ef(ji))/(4*k*T)))))/(e^2/h_plank); % real part of conductivity
        cond_gi(ji,:)=((4*Ef(ji))* (e^2/(4*h_plank))./(h_plank*omega*pi)*(1-2*Ef(ji)^2/(9*t^2))-(1+(h_plank*omega).^2/(36*t^2))* e^2/(4*h_plank)/pi.*log(abs(h_plank*omega+2*Ef(ji))./abs(h_plank*omega-2*Ef(ji))))/(e^2/h_plank); % imaginary part of conductivity
        cond_g(ji,:)=cond_gr(ji,:)+cond_gi(ji,:)*1i; %e^2/h_plank
        %perm(ji,:)=(5.5+(1i*cond_g(ji,:)*(e^2/h_plank)./(dgr*omega))*7.05*10^29); % permittivity, 7.05*10^29 is the convertion number from e^2/h_plank to S over the permittivity of vacuum
        perm(ji,:)=(5.5/(dgr/(0.34*10^-9))+(1i*cond_g(ji,:)*(e^2/h_plank)./(dgr*omega))*7.05*10^29); % permittivity, 7.05*10^29 is the convertion number from e^2/h_plank to S over the permittivity of vacuum
    else 
        cond_gi(ji,:)=( 1i*e^2/(4*pi*h_plank)  *  log(  (2*abs(Ef(ji))-h_plank*omega-1i*Gama)  ./  (2*abs(Ef(ji))+h_plank*omega+1i*Gama)  ))/(e^2/h_plank);
        cond_gr(ji,:)=( 1i*e^2./(pi*h_plank*(h_plank*omega+1i*Gama))  .*  (Ef(ji)+2*k*T*log(exp(-Ef(ji)/(k*T))+1)))/(e^2/h_plank);
        cond_g(ji,:)=(cond_gr(ji,:)+cond_gi(ji,:));  %e^2/h_plank
        %perm(ji,:)=(5.5+(1i*cond_g(ji,:)*(e^2/h_plank)./(dgr*omega))*7.05*10^29); % permittivity, 7.05*10^29 is the convertion number from e^2/h_plank to S over the permittivity of vacuum
        perm(ji,:)=(5.5/(dgr/(0.34*10^-9))+(1i*cond_g(ji,:)*(e^2/h_plank)./(dgr*omega))*7.05*10^29); % permittivity, 7.05*10^29 is the convertion number from e^2/h_plank to S over the permittivity of vacuum
    end
    %perm(ji,:)=(5.5+(1i*cond_g(ji,:)*(e^2/h_plank)./(dgr*omega))*7.05*10^29); % permittivity, 7.05*10^29 is the convertion number from e^2/h_plank to S over the permittivity of vacuum
    permlorenrz(ji,:)=FL*EL^2./((EL^2-(h_plank*omega).^2).^2+GL^2.*(h_plank*omega).^2).*((EL^2-(h_plank*omega).^2)+1i*GL.*h_plank*omega); % Lorentz term
    ng(ji,:)=sqrt(perm(ji,:)); % refractive index    
end
%--------------------
nr=real(ng);
ni=imag(ng);
%permr=real(perm);
%permi=imag(perm);
%condr=real(cond_g);
%condi=imag(cond_g);           
% *********************************************************
%Export data to text
% *********************************************************
Graphene=table(ene2',nr',ni');
writetable(Graphene,'MNPBEM17/Material/@epstable/Graphene.dat','Delimiter',' ','WriteVariableNames',0)
%output messages-------------
%sprintf('The Fermi level is %0.2f \n',Ef) 
%sprintf('The temperature is %0.2f \n',T)
%sprintf('The minimum wavelength is %0.2f \n',min(lambda*10^9))
%sprintf('The maximum wavelength is %0.2f \n',max(lambda*10^9))
%sprintf('The number of wavelength steps is %0.1f \n',length(nr))  
end



