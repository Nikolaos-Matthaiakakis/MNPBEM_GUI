function [p, op] = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
% Al Heptamer
% Eps table
epstab = { epsconst( 1 ), epstable( 'gold.dat' )};%  table of dielectric functions

% Particle
len = [ d_length, d_length ]; %  Dimensions of particle
poly = polygon( 20, 'size', len, 'dir', -1 );%  polygon for disk
poly_sub = polygon( 4, 'size', [ d_length3, d_length3  ], 'dir', 1 ); % substrate
edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
hdata = struct( 'hmax', hdata_tmp ); %  hdata
h_vert=(1/2) * sqrt(3) * shift_l;% vertical distance
x_vert=shift_l/2;% longitudinal distance
p1 = shift( poly, [0 , 0] );
p2 = shift( poly, [-shift_l , 0] );        
p3 = shift( poly, [shift_l , 0 ] );
p4 = shift( poly, [x_vert , h_vert ] );
p5 = shift( poly, [-x_vert , -h_vert ] );
p6 = shift( poly, [x_vert , -h_vert ] );
p7 = shift( poly, [-x_vert , h_vert ] );
p = tripolygon( [poly_sub,p1,p2,p3,p4,p5,p6,p7], edge, 'hdata', hdata); %  Disk ,p1,p2,p3,p4,p5,p6,p7
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
% Rotate
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
end





