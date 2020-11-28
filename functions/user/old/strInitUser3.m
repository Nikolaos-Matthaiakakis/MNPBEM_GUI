%% Selectrion of nanostructure
function p = strInit(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en)
typeStr=5;
if typeStr==0 %  Nanosphere
    if substrate==0
        p=trisphere( grid_1, d_length, 'interp', 'curv' );
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
        p=trisphere( grid_1, d_length, 'interp', 'curv' );
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - d_length/2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------%    
elseif typeStr==1 %  Nanodisk  
    if substrate==0
        poly = polygon( 25, 'size', [ 1, 1 ] * d_length ); %  Polygon for disk
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p= tripolygon( poly, edge, 'hdata', hdata );
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
        poly = polygon( 25, 'size', [ 1, 1 ] * d_length ); %  Polygon for disk
        edge = edgeprofile( edge_profile, grid_2, 'mode', '01', 'min', 0  ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        [p,poly]= tripolygon( poly, edge, 'hdata', hdata );            
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
%-------------------------%    
elseif typeStr==2 %  Nanotriangle
    if substrate==0
        poly = round( polygon( 3, 'size', [ d_length,  2 / sqrt( 3 ) * d_length ] ) ); %  Polygon for triangle
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for triangle
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( poly, edge, 'hdata', hdata); %  Extrude polygon to nanoparticle            
        if Double_en==1
            p1t = rot( p, 180 , [0,0,1] ); %  Rotation around z (degrees)
            p1 = shift( p1t, [ -d_length/2-shift_l , 0, 0 ] );            
            p2 = shift( p, [ d_length/2+shift_l, 0, 0 ] );            
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
                if qd_en==1 && Source_op==3
                    p2t=trisphere( grid_1, 15, 'interp', 'curv' );
                    p2 = shift( p2t, [ dip_p(1), dip_p(2), dip_p(3)] );
                    p = comparticle( epstab, { p, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
                else
                    p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
                end        
            end
        end
    else
        poly = round( polygon( 3, 'size',[ d_length,  2 / sqrt( 3 ) * d_length ] )); %  Polygon for disk
        edge = edgeprofile( edge_profile, grid_2, 'mode', '01', 'min', 0  ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        [p,poly]= tripolygon( poly, edge, 'hdata', hdata );            
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
%-------------------------%    
elseif typeStr==3 %  Polygon
    if substrate==0
        poly = round( polygon( angles_n, 'size', [ d_length,  2 / sqrt( 3 ) * d_length ] ) ); %  Polygon
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for polygon
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( poly, edge, 'hdata', hdata ); %  Extrude polygon to nanoparticle
        p = shift( p, [ 15, 0, 0 ] );
        if Double_en==1
            p1t = rot( p, 180 , [0,0,1] ); %  Rotation around z (degrees)
            p1 = shift( p1t, [ -d_length/2-shift_l/2 , 0, 0 ] );            
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
         end  
    else
        poly = round( polygon( angles_n, 'size',[ d_length,  2 / sqrt( 3 ) * d_length ] )); %  Polygon for disk
        edge = edgeprofile( edge_profile, grid_2, 'mode', '01', 'min', 0  ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        [p,poly]= tripolygon( poly, edge, 'hdata', hdata );            
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
%-------------------------%    
elseif typeStr==4 %  Nanorod vertical    
    if substrate==0
        p = trirod( edge_profile, d_length, [ grid_2, grid_2, grid_1 ] ); %  Initialize nanorod
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
        p = trirod( edge_profile, d_length, [ grid_2, grid_2, grid_1 ] ); %  Initialize nanorod
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - d_length/2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------%    
elseif typeStr==5 %  Nanorod horizontal
    if substrate==0
        len = [ d_length, d_length2 ]; %  Dimensions of particle
        poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( poly, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon
        len = [ d_length3, d_length3 ]; %  Dimensions of particle
        poly = round( polygon( 4, 'size', len ), 'rad', 5 ); %  Polygon
        edge = edgeprofile( sub_t, grid_2 ); %  Edge profile
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p2 = tripolygon( poly, edge, 'hdata', hdata, 'interp', 'curv' ); %  Extrude polygon
        if Double_en==1
            p1 = p;
            p2 = shift( p2, [ 0 , 0, -sub_t/2-edge_profile/2-0.5 ] );
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
elseif typeStr==6 %  Nano ellipsoid
    if substrate==0            
        ax = [ d_length, d_length2 ,edge_profile];
        p = scale( trisphere( grid_1, 1 , 'interp', 'curv' ), ax );
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
        ax = [ d_length, d_length2 ,edge_profile];
        p = scale( trisphere( grid_1, 1 ), ax );
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - d_length*2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------%    
elseif typeStr==7 %  Cube
    if substrate==0
        p=tricube( grid_1, d_length ); 
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
        p=tricube( grid_1, d_length );
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - d_length/2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------% 
elseif typeStr==8 %  Torus  
    if substrate==0
        poly1 = polygon( 30, 'size', [ d_length, d_length ], 'dir', 1 );
        poly2 = polygon( 30, 'size', [ d_length2, d_length2 ], 'dir',  -1 );
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( [ poly1, poly2 ], edge, 'hdata', hdata  );%  extrude particle 
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
        poly1 = polygon( 30, 'size', [ d_length, d_length ], 'dir', 1 );
        poly2 = polygon( 30, 'size', [ d_length2, d_length2 ], 'dir',  -1 );
        edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for disk
        hdata = struct( 'hmax', hdata_tmp ); %  hdata
        p = tripolygon( [ poly1, poly2 ], edge, 'hdata', hdata  );%  extrude particle
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - edge_profile/2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------%    
elseif typeStr==9 %  Torus 2   
    if substrate==0
        p = tritorus( d_length, edge_profile, [ grid_1, grid_2 ] ); 
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
        p = tritorus( d_length, edge_profile, [ grid_1, grid_2 ] ); 
        x = d_length3 * linspace( - 0.5, 0.5, 21 );
        [ verts, faces ] = fvgrid( x, x );            
        up = flipfaces( shift( particle( verts, faces ), [ 0, 0, - edge_profile/2 ] ) );%  make particle
        lo = shift( up, [ 0, 0, - sub_t ] );
        %up=flipfaces(up);
        if substrate2==0
            p = comparticle( epstab, { p, up, lo },[ 2, 1; 1, 3; 3, 1 ], [ 1, 2 ], op );%  make particle
        else
            lo2 = shift( lo, [ 0, 0, - sub_t2] );
            p = comparticle( epstab, { p, up, lo, lo2 },[ 2, 1; 1, 3; 3, 4; 4, 1 ], [ 1, 2 ], op );%  make particle
        end
        p = flip( p, 3 );
    end
%-------------------------% 
elseif typeStr==10 %  Nanohole
    poly1 = polygon( 4, 'size', [ d_length2, d_length2  ], 'dir', 1 );
    poly2 = round(polygon( angles_n, 'size', [ d_length, d_length ], 'dir',  -1 ));
    edge = edgeprofile( edge_profile, grid_2 ); %  Edge profile for disk
    hdata = struct( 'hmax', hdata_tmp ); %  hdata
    p = tripolygon( [ poly1, poly2 ], edge, 'hdata', hdata  );%  extrude particle 
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
end
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
figure
plot(p)
end