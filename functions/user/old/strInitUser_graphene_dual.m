%% Graphene dual structure
function p = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
        %shift_l=15;
        len = [ d_length, d_length2 ]; %  Dimensions of particle
        len_red = [ d_length/1.5, d_length2/1.5 ]; %  Dimensions of particle
        poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
        poly_red = round( polygon( 4, 'size', len_red ), 'rad', 5 ); %  Polygon
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( poly, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon
        p_red = tripolygon( poly_red, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon
            p1 = shift( p, [ -d_length/2-shift_l/2 , -d_length2/2-shift_l/2, 0 ] );          
            ptemp = rot( p_red, 90 , [0,0,1] ); %  Rotation around z (degrees)
            p2 = shift( ptemp, [ (len_red(2))/2+shift_l/2, len_red(1)/2+shift_l/2, 0 ] );            
            p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
        p = flip( p, 3 );
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
end