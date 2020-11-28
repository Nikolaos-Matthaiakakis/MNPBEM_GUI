function [p, op] = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
% Al Heptamer
% Eps table
if Cover_en==1 % currently only for single disk
    epstab = { epsconst( 1 ), epstable( 'Aluminum.dat' ), epsconst( 3.168 ) };% table of dielectric functions Al2O3
else
    epstab = { epsconst( 1 ), epstable( 'Aluminum.dat' )};%  table of dielectric functions
end
% Particle
len = [ d_length, d_length ]; %  Dimensions of particle
poly = polygon( 20, 'size', len );%  polygon for disk
edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
hdata = struct( 'hmax', hdata_tmp ); %  hdata
h_vert=(1/2) * sqrt(3) * shift_l;% vertical distance
x_vert=shift_l/2;% longitudinal distance
p = tripolygon( poly, edge, 'hdata', hdata); %  Disk
%p=trisphere( grid_1, d_length, 'interp', 'curv' ); %Sphere
if Cover_en==1 % currently only for single disk
    p1 = shift( p, [0 , 0, 0 ] );     
    p_cover1 = coverlayer.shift( p1, d_c );      
    p2 = shift( p, [-shift_l , 0, 0 ] );  
    p_cover2 = coverlayer.shift( p2, d_c );
    p3 = shift( p, [shift_l , 0, 0 ] );
    p_cover3 = coverlayer.shift( p3, d_c );
    p4 = shift( p, [x_vert , h_vert, 0 ] );
    p_cover4 = coverlayer.shift( p4, d_c );
    p5 = shift( p, [-x_vert , -h_vert, 0 ] );
    p_cover5 = coverlayer.shift( p5, d_c );
    p6 = shift( p, [x_vert , -h_vert, 0 ] );
    p_cover6 = coverlayer.shift( p6, d_c );
    p7 = shift( p, [-x_vert , h_vert, 0 ] );  
    p_cover7 = coverlayer.shift( p7, d_c );
    p = comparticle( epstab, { p1, p2, p3, p4, p5, p6, p7, p_cover1, p_cover2, p_cover3, p_cover4, p_cover5, p_cover6, p_cover7 }, [ 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; ], 1, 2, 3, 4, 5, 6, 7, op );
else
    p1 = shift( p, [0 , 0, 0 ] );          
    p2 = shift( p, [-shift_l , 0, 0 ] );       
    p3 = shift( p, [shift_l , 0, 0 ] );
    p4 = shift( p, [x_vert , h_vert, 0 ] );
    p5 = shift( p, [-x_vert , -h_vert, 0 ] );
    p6 = shift( p, [x_vert , -h_vert, 0 ] );
    p7 = shift( p, [-x_vert , h_vert, 0 ] );    
    p = comparticle( epstab, { p1, p2, p3, p4, p5, p6, p7 }, [ 2, 1; 2, 1; 2, 1; 2, 1; 2, 1; 2, 1; 2, 1; ], 1, 2, 3, 4, 5, 6, 7, op );
end
% Rotate
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
end





