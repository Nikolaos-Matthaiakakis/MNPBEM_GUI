function [ff]=farfield_qd(op, enei, ien, dip_d, spec)
%%  initialization
%  table of dielectric functions
epstab = { epsconst( 1 ), epstable( 'lorentz.dat' ) };
%  diameter of sphere
diameter = 15;
%  initialize sphere
p = comparticle( epstab, { trisphere( 144, diameter ) }, [ 2, 1 ], 1, op );

%%  dipole oscillator
%  positions of dipole
z = 0;
%  compoint
pt = compoint( p, [ 0 * z, 0 * z, z ] );
%  dipole excitation
if dip_d==1
    dip  = dipole( pt, [ 1, 0, 0], op ); %  Dipole excitation X
elseif dip_d==2
    dip  = dipole( pt, [ 0, 1, 0], op ); %  Dipole excitation Y
else
    dip  = dipole( pt, [ 0, 0, 1], op ); %  Dipole excitation Z
end   

%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei(ien) );
%  Farfield
ff=farfield(spec,sig);

end