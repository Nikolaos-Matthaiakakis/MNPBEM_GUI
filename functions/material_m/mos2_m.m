% *********************************************************
% MoS2 optical model (d=0.65nm)
% *********************************************************
%Bablu Mukherjee, Frank Tseng, Daniel Gunlycke, 
%Kiran Kumar Amara, Goki Eda, and Ergun Simsek, 
%"Complex electrical permittivity of the monolayer 
%molybdenum disulfide (MoS2) in near UV and visible,
%" Opt. Mater. Express 5, 447-455 (2015) 

function mos2_m( )
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
Vg=-10; %  V
fermi=0; %  1 enables calculation from voltage

%----------------------------------------------------------
%Main program**********************************************
%----------------------------------------------------------

% *********************************************************
% Constants and variables
% *********************************************************
c = 2.99792458*10^8; %  Speed of light (m/s)
h_p=6.582*10^-16; %  eV.s
h_p2=1.0546*10^-34; %  J.s 
e=1.6*10^-19; %  Electron charge (C)
Kb=8.6*10^-5; %  Boltzmann constant (eVK^-1 )
% *********************************************************
% Oscillator settings
% *********************************************************
j=6; %  Number of LD oscillators
aj=[2.0089*10^5 5.7534*10^4 8.1496*10^4 8.2293*10^4 3.3130*10^5 4.3906*10^6]/40; %  Oscillator strength 
bj=[1.0853*10^-2 5.9099*10^-2 1.1302*10^-1 1.1957*10^-1 2.8322*10^-1 7.8515*10^-1]; %  Damping coefficient (eV)
E_wj=[0 1.88 2.03 2.78 2.91 4.31]; %  Resonance frequency (Hz)
E_wp=28.3/1000; %  Plasma frequency (eV)
a=23.224; %  Maximum G
u=2.7723; %  mean G
sigma=0.3089; %  variance
eps_inf=4.44; %  Background permittivity
% *********************************************************
% Wavelength
% *********************************************************
ene2= linspace( minene2, maxene2, stepene2 ); %  eV
%ene2=1.5:0.01:3.5; 
if length(ene2)==1
    ene2=0.4:0.01:ene2+2; %  eV
end
% *********************************************************
% Energy Conversion
% *********************************************************
lambda=((4.136*10^-15)*3*10^8)./ene2; %m
omega=2*pi*c./lambda; %hz
lambdap=((4.136*10^-15)*3*10^8)./E_wp; %m
omegap=2*pi*c./lambdap; %hz
% *********************************************************
% Fermi level
% *********************************************************
Ef=Ef_1;
if fermi==1
    Vos=Vg+107; %  V
    C=1.2*10^-8;% F*m^-2
    ne=abs(C*Vos/e*10^4); %  m^-2
    me=0.35*9.10938356*10^-31; %  kg 0.511%
    Ef=h_p2^2*pi*ne/(2*me*e^2)*e; %  Fermi level eV
end
% *********************************************************
% Refractive index of MoS2
% *********************************************************
delta=e^-(12*pi*(Ef-Kb*T).^2); %  Empirical term for gate
for n=1:6
    %eps_LD(n,:)=eps_inf+aj(n)/delta*omegap^2./(wj(n)^2-omega.^2-1i*omega.*bj(n)); %  DL model
    eps_LD(n,:)=aj(n)/delta*E_wp^2./((E_wj(n))^2-ene2.^2-1i*ene2.*(bj(n))*delta); %  DL model
end
eps_LD_tot=eps_inf+eps_LD(1,:)+eps_LD(2,:)+eps_LD(3,:)+eps_LD(4,:)+eps_LD(5,:)+eps_LD(6,:); %  LD sum permittivity
eps_gi=a*exp(-(h_p*omega-u).^2/(2*sigma^2)); %  G model imag
eps_gr=kkrebook2(ene2,eps_gi,0); %  K-K for real part of G model
eps=eps_LD_tot+eps_gr+1i*eps_gi; %  Final permittivity
ng=sqrt(eps); %  Refractive index  
%--------------------
nr=real(ng); %  n
ni=imag(ng); %  k
% *********************************************************
%Export data to text
% *********************************************************
mos2=table(ene2',nr',ni');
writetable(mos2,'MNPBEM17/Material/@epstable/mos2.dat','Delimiter',' ','WriteVariableNames',0)

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
%figure
%plot(lambda*10^9, imag(eps_LD(1,:)),lambda*10^9, imag(eps_LD(2,:)),lambda*10^9, imag(eps_LD(3,:)),lambda*10^9, imag(eps_LD(4,:)),lambda*10^9, imag(eps_LD(5,:)),lambda*10^9,imag( eps_LD(6,:)),lambda*10^9,eps_gi,lambda*10^9,imag(eps_LD_tot)+eps_gi)
%figure
%plot(lambda*10^9, real(eps_LD(1,:)),lambda*10^9, real(eps_LD(2,:)),lambda*10^9, real(eps_LD(3,:)),lambda*10^9, real(eps_LD(4,:)),lambda*10^9, real(eps_LD(5,:)),lambda*10^9,real( eps_LD(6,:)),lambda*10^9,eps_gr,lambda*10^9,real(eps_LD_tot)+eps_gr)
%figure
%plot(lambda*10^9, imag(eps),lambda*10^9,real(eps))
end
