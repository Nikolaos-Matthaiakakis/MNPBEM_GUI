%% User defined nanostructure
function p = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
    %  p returns the structure for simulation
    %  Input values can be used as in function strInit to allow for control
    %  through the GUI
    %
    %  Plot areas are still defined through the GUI, thus if the GUI input
    %  values here are not used, defining the structure size in the GUI is
    %  still required!!!     
    if substrate==0
        len = [ d_length, d_length2 ]; %  Dimensions of particle
        poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( poly, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon
        if Double_en==1
            p1 = shift( p, [ -d_length/2-shift_l/2 , 0, 0 ] );
            p2 = shift( p, [ d_length/2+shift_l/2, 0, 0 ] );
            if qd_en==1 && Source_op==3
                p3t=trisphere( grid_1, 15, 'interp', 'curv' );
                p3 = shift( p3t, [ dip_p(1), dip_p(2), dip_p(3)] );
                p = comparticle( epstab, { p1, p2, p3 }, [ 2, 1; 3, 1; 4, 1 ], 1, 2, op );
            else
                p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
            end
        elseif Cover_en==1
            p1 = p;
            p2 = coverlayer.shift( p1, d_c );
            if qd_en==1 && Source_op==3
                p3t=trisphere( grid_1, 15, 'interp', 'curv' );
                p3 = shift( p3t, [ dip_p(1), dip_p(2), dip_p(3)] );
                p = comparticle( epstab, { p1, p2, p3 }, [ 2, 3; 3, 1; 4, 1 ], 1, 2, op );
            else
                p = comparticle( epstab, { p1, p2 }, [ 2, 3; 3, 1 ], 1, 2, op ); 
            end                       
        else                
            if qd_en==1 && Source_op==3
                p2t=trisphere( grid_1, 15, 'interp', 'curv' );
                p2 = shift( p2t, [ dip_p(1), dip_p(2), dip_p(3)] );
                p = comparticle( epstab, { p, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
            else
                p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
            end
        end
    else
        len = [ d_length, d_length2 ]; %  Dimensions of particle
        poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
        edge = edgeprofile( edge_profile, grid_2, 'mode', '01', 'min', 0 ); %  Edge profile
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        [p,poly] = tripolygon( poly, edge, 'hdata', hdata ); %  Extrude polygon           
        %----------- Create substrate
        [ pup, plo ] = select( p, 'carfun', @( x, y, z ) z > 1e-3 );%  split init upper and lower part            
        poly2 = polygon3( polygon( 4, 'size', [ d_length3, d_length3 ] ), 0 );%  polygon for plate            
        up = plate( [ poly2, set( poly, 'z', 0 ) ], 'dir', - 1 );%  upper plate
        %  lower plate, we use FVGRID to get quadrilateral face elements
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        lo = flipfaces( shift( particle( verts, faces ), [ 0, 0, - sub_t ] ) );%  make particle
        if substrate2==0
            p = comparticle( epstab, { pup, plo, up, lo },[ 2, 1; 2, 3; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { pup, plo, up, lo, lo2 },[ 2, 1; 2, 3; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
end