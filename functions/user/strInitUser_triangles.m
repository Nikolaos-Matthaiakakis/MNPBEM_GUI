function [p, op] = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
% Al Heptamer
% Eps table
%epstab = { epsconst( 1 ), epstable( 'gold.dat' )};%  table of dielectric functions
epstab = { epsconst( 1 ), epstable( 'gold.dat' ),epstable( 'mos2.dat' )};%  table of dielectric functions
% Particle
len = [ d_length, d_length* 2 / sqrt( 3 ) ]; %  Dimensions of particle
poly = round(polygon( 3, 'size', len ));%  polygon for disk
poly_sub = round(polygon( 4, 'size', d_length3 ));%  polygon for disk
edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
edge_sub = edgeprofile( sub_t, grid_2 ); %  Edge profile
hdata = struct( 'hmax', hdata_tmp ); %  hdata
h_vert=(1/2) * sqrt(3) * shift_l;% vertical distance
x_vert=shift_l/2;% longitudinal distance
p = tripolygon( poly, edge, 'hdata', hdata); %  Disk        
p2 = shift( p, [-shift_l , 0, 0 ] );
p3=rot(p,180);
p3 = shift( p3, [shift_l , 0, 0 ] );
p4 = shift( p, [x_vert , h_vert, 0 ] );
p5=rot(p,180);
p5 = shift( p5, [-x_vert , -h_vert, 0 ] );
p6 = shift( p, [x_vert , -h_vert, 0 ] );
p7=rot(p,180);
p7 = shift( p7, [-x_vert , h_vert, 0 ] ); 
p_sub = tripolygon( poly_sub, edge_sub, 'hdata', hdata); %  Disk
p_sub=shift(p_sub, [0 , 0, -edge_profile/2-sub_t/2 ]);
%p = comparticle( epstab, {  p2, p3, p4, p5, p6, p7 }, [ 2, 1; 2, 1; 2, 1; 2, 1; 2, 1; 2, 1 ], 1, 2, 3, 4, 5, 6, op );
p = comparticle( epstab, {  p2, p3, p4, p5, p6, p7, p_sub }, [ 2, 1; 2, 1; 2, 1; 2, 1; 2, 1; 2, 1 ; 3, 1], 1, 2, 3, 4, 5, 6, 7, op );
% Rotate
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
end