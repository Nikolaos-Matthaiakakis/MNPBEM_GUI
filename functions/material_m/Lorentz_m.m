% *********************************************************
% MoS2 optical model (d=0.65nm)
% *********************************************************
%Bablu Mukherjee, Frank Tseng, Daniel Gunlycke, 
%Kiran Kumar Amara, Goki Eda, and Ergun Simsek, 
%"Complex electrical permittivity of the monolayer 
%molybdenum disulfide (MoS2) in near UV and visible,
%" Opt. Mater. Express 5, 447-455 (2015) 

function Lorentz_m( )
% *********************************************************
% Code inputs
% *********************************************************
% Load input variables  
mat = dir('tmp/*.mat');
for q = 1:length(mat) 
    load(sprintf(['tmp/' mat(q).name])); 
end 

%----------------------------------------------------------
%Main program**********************************************
%----------------------------------------------------------

% *********************************************************
% Constants and variables
% *********************************************************
c = 2.99792458*10^8; %  Speed of light (m/s)
bj=gama_eg; %  Damping coefficient (eV)
E_wj=omega_eg; %  Resonance frequency (Hz)
E_wp=omega_pl; %  Plasma frequency (eV)
ene2= linspace( minene2, maxene2, stepene2 ); %  eV
if length(ene2)==1
    ene2=0.4:0.01:ene2+2; %  eV
end

% *********************************************************
% Refractive index
% *********************************************************
eps_LD=E_wp^2./((E_wj)^2-ene2.^2-1i*ene2.*(bj)); %  L model
ng=sqrt(1+eps_LD(1,:)); %  Refractive index  
%--------------------
nr=real(ng); %  n
ni=imag(ng); %  k
% *********************************************************
%Export data to text
% *********************************************************
lorentz=table(ene2',nr',ni');
writetable(lorentz,'MNPBEM17/Material/@epstable/lorentz.dat','Delimiter',' ','WriteVariableNames',0)

% *********************************************************
%Test
% *********************************************************
%eps_final=kkrebook2(ene2,eps,0);
%figure
%plot(ene2, imag(eps),ene2, imag(eps_LD),'--',ene2, eps_gi,'.')
%figure
%plot(ene2, real(eps),ene2, real(eps_LD),'--',ene2, eps_gr,'.')
%figure
%plot(lambda*10^9, imag(eps),lambda*10^9, imag(eps_LD),'--',lambda*10^9, eps_gi,'.')
%figure
%plot(lambda*10^9, real(eps),lambda*10^9, real(eps_LD),'--',lambda*10^9, eps_gr,'.')
%figure
%plot(lambda*10^9, nr,lambda*10^9, ni,'--')
end