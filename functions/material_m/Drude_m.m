% *********************************************************
% ITO optical model 
% *********************************************************
%Matthaiakakis, Nikolaos      (2017)     
%Dynamic modulation of plasmon excitations in monolayer graphene.  
%University of Southampton, Doctoral Thesis, 231pp. 
function Drude_m( )
% *********************************************************
% Code inputs
% *********************************************************
% Load input variables  
mat = dir('tmp/*.mat');
for q = 1:length(mat) 
    load(sprintf(['tmp/' mat(q).name])); 
end 
%T=300; % Temperature K

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
% *********************************************************
% Wavelength convert
% *********************************************************
lambda_optmp=(h_plankeV*c)./omega_p; %um
lambda_op=lambda_optmp*10^-6; %m
omegap=2*pi*c./lambda_op;
%-----
ene2= linspace( minene2, maxene2, stepene2 ); %  eV
if length(ene2)==1
    ene2=[ene2 ene2+0.01];
end
la= (h_plankeV*c)./ene2; %wavelength range in um
lambda=la*10^-6; %m
omega_=2*pi*c./lambda;
%-----
gama_=(h_plankeV*c)/gamaeV;
lambda_G=gama_*10^-6; %m
gama=2*pi*c./lambda_G;
% *********************************************************
% Permittivity and Refractive index
% *********************************************************
perm=eps_inf-omegap.^2./(omega_.^2+1i.*gama.*omega_);
% *********************************************************
%Export data to text
% *********************************************************
n=sqrt(perm);
nr=real(n);
ni=imag(n);
Drude_nk=table(ene2',nr',ni');
writetable(Drude_nk,'MNPBEM17/Material/@epstable/Drude.dat','Delimiter',' ','WriteVariableNames',0)
% *********************************************************
%Test
% *********************************************************
%figure
%plot(ene2,nr,ene2,ni,'--')