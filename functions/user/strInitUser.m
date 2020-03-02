%% Au Graphene hybrid structure
function p = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
% Gold rectangle
len = [ d_length, d_length2 ]; %  Dimensions of particle
poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
hdata = struct( 'hmax', hdata_tmp ); %  hdata
p1 = tripolygon( poly, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon             
% Graphene disk
poly2 = polygon( 25, 'size', [ 1, 1 ] * d_length/6 ); %  Polygon for disk
edge2 = edgeprofile( 1, grid_2 ); %  Edge profile for disk
hdata2 = struct( 'hmax', hdata_tmp ); %  hdata
p2temp= tripolygon( poly2, edge2, 'hdata', hdata2 );
p2 = shift( p2temp, [d_length/2+d_length/12,d_length2/2+d_length/12, 0 ] ); 
% Gold rectangle
p3=shift( p1, [d_length+d_length/6,d_length2+d_length/6, 0 ] );
% Final structure        
p = comparticle( epstab, { p1, p2, p3 }, [ 2, 1; 3, 1; 2, 1 ], 1, 2, 3, op );
%p = comparticle( epstab, { p2}, [ 2, 1 ], 1, op );
% Rotation
p = flip( p, 3 );
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)