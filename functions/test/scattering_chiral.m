function [ sca, dsca ] = scattering_chiral( field1, field2, cp, medium )
%  SCATTERING - Radiated power for electromagnetic fields.
%
%  Usage  :
%    [ sca, dsca ] = scattering( field )
%  Input
%    field      :  compstruct object with scattered electromagnetic fields
%    medium     :  compute total radiated pwoer only in given medium
%  Output
%    sca        :  total radiated power
%    dsca       :  differential radiated power for particle surface
%                  given in FIELD

%  particle surface for fields at 'infinity'
pinfty = field1.p;
%  scattered electric and magnetic fields
ep = field1.e; 
hp = field1.h; 
es = field2.e; 
hs = field2.h; 
%  Poynting vector in direction of outer surface normal
%dsca = 0.5 * real( inner( pinfty.nvec, cross( ep, conj( hp ), 2 ) ) );
%size(dsca)
if cp==1
    dsca = 0.25 * real( inner( pinfty.nvec, abs((conj(es)-1i*conj(ep)).^2 )) );
else
    dsca = 0.25 * real( inner( pinfty.nvec,  cross((conj(es)+1i*conj(ep)),(conj(es)+1i*conj(ep)),2) ) );
end
%size(dsca1)
%  area of boundary elements
area = pinfty.area;
%  total cross section only in given medium
if exist( 'medium', 'var' ) && ~isempty( medium )
  area( pinfty.expand( num2cell( pinfty.inout( :, end ) ) ) == medium ) = 0;
end
%  total radiated power
sca = squeeze( matmul( reshape( area, 1, [] ), dsca ) );
%  convert differential radiated power to compstruct object
dsca = compstruct( pinfty, field1.enei, 'dsca', dsca );