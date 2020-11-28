function [p, op, exc_flag, bem] = strInitUser_Al_comp(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata_tmp, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en, ene2, Field_en, xe, ye)
%% Al Heptamer with substrate
units;
%% enable
disk_set=1;% switch between disk and heptamer
field_set=1;% enable field calculation
%% input
angle_set=0;% angle of incident light
%epstab = { epsconst( 1 ), epstable( 'Drude.dat' ), epsconst( 3.168 ), epsconst( 3.168 ), epstable( 'Aluminum.dat' )};% table of dielectric functions
epstab = { epsconst( 1 ), epstable( 'Drude.dat' ), epsconst( 3.168 ), epsconst( 3.168 )};% table of dielectric functions
[ x, z ] = meshgrid( linspace( - 150, 150, 100 ), linspace( - 80, 180, 100 ) );% Field area
%% Init
enei = eV2nm ./ ene2;
op = layerstructure.options; %  default options for layer structure
ztab = [ sub_t2, 0 ];        
%layer = layerstructure( epstab, [ 1, 4, 5 ], ztab, op );%  set up layer structure
layer = layerstructure( epstab, [ 1, 3, 4 ], ztab, op );%  set up layer structure
op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'layer', layer ); % Options for BEM (Retarded)
%% Particle
len = [ d_length, d_length ]; %  Dimensions of particle
poly = polygon( 20, 'size', len );%  polygon for disk
edge = edgeprofile( edge_profile, grid_2,'mode', '01', 'min', 1e-3 ); %  Edge profile
%    MODE '01' produces a rounded edge on top and a sharp edge on bottom,
%    MIN controls the lower z-value of the nanoparticle 
hdata = struct( 'hmax', hdata_tmp ); %  hdata
h_vert=(1/2) * sqrt(3) * shift_l;% vertical distance
x_vert=shift_l/2;% longitudinal distance
p = tripolygon( poly, edge, 'hdata', hdata); %  Disk
%p=trisphere( grid_1, d_length, 'interp', 'curv' ); %Sphere
p_c = scale(p,2);
%p_c = scale(p,(d_length+d_c)/d_length);
p = shift( p, [ 0, 0, sub_t2] );%  shift nanosphere above substrate
p_c = shift( p_c, [ 0, 0, sub_t2] );%  shift nanosphere above substrate
if Cover_en==1 % currently only for single disk
    if disk_set==1
        p1 = shift( p, [0 , 0, d_c+0.1+20 ] );
        p_c1= shift( p_c, [0 , 0, d_c] );
        p = comparticle( epstab, { p1, p_c1 }, [ 2, 3; 3, 1], 1, 2, op );% middle disk only    
    else
        p1 = shift( p, [0 , 0, d_c ] );     
        p_c1 = shift( p_c, [0 , 0, d_c ] );     
        p2 = shift( p, [-shift_l , 0, d_c ] );  
        p_c2 = shift( p_c, [-shift_l , 0, d_c ] );
        p3 = shift( p, [shift_l , 0, d_c ] );
        p_c3 = shift( p_c, [shift_l , 0, d_c ] );
        p4 = shift( p, [x_vert , h_vert, d_c ] );
        p_c4 = shift( p_c, [x_vert , h_vert, d_c ] );
        p5 = shift( p, [-x_vert , -h_vert, d_c ] );
        p_c5 = shift( p_c, [-x_vert , -h_vert, d_c ] );
        p6 = shift( p, [x_vert , -h_vert, d_c ] );
        p_c6 = shift( p_c, [x_vert , -h_vert, d_c ] );
        p7 = shift( p, [-x_vert , h_vert, d_c ] );  
        p_c7 = shift( p_c, [-x_vert , h_vert, d_c ] );
        p = comparticle( epstab, { p1, p2, p3, p4, p5, p6, p7, p_c1, p_c2, p_c3, p_c4, p_c5, p_c6, p_c7 }, [ 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 2, 3; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; 3, 1; ], 1, 2, 3, 4, 5, 6, 7, op );
    end
else
    if disk_set==1
        p1 = shift( p, [0 , 0, 0 ] );  
        p = comparticle( epstab, { p1}, [ 2, 1], 1, op );% middle disk only  
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
end
%% Rotate
p = rot( p, xtilt , [1,0,0] ); %  Rotation around x (degrees)
p = rot( p, ytilt , [0,1,0] ); %  Rotation around y (degrees)
p = rot( p, ztilt , [0,0,1] ); %  Rotation around z (degrees)
%%  light propagation angles
%theta = pi / 180 * reshape( linspace( 0, 80, 5 ), [], 1 );
theta = pi / 180 * reshape( linspace( angle_set, angle_set, 1 ), [], 1 );
%  TM mode, excitation from above
dir = [ sin( theta ), 0 * theta, - cos( theta ) ];
pol = [ cos( theta ), 0 * theta,   sin( theta ) ];
%% Green function
%  grid where fields will be computed
%[ x, z ] = meshgrid( linspace( - 30, 30, 81 ), linspace( - 30, 60, 101 ) );
%  make compoint object
%    it is important that COMPOINT receives the OP structure because it has
%    to group the points within the layer structure
if field_set==1
    pt = compoint( p, [ x( : ),  0 *z( : ), z( : ) ], op );
    %  Table creation
    if ~exist( 'greentab', 'var' ) || ~greentab.ismember( layer, enei, p, pt )
        %  automatic grid for tabulation
        %    we use a rather small number NZ for tabulation to speed up the
        %    simulations
        tab = tabspace( layer, p, pt );
        %  Green function table
        greentab = compgreentablayer( layer, tab );
        %  precompute Green function table
        %    for a more accurate simulation of the layer the number of
        %    wavelenghts should be increased
        greentab = parset( greentab, linspace( min(enei), max(enei), length(enei) ), op );
    end  
else
    if ~exist( 'greentab', 'var' ) || ~greentab.ismember( layer, enei, p )
        %  automatic grid for tabulation
        %    we use a rather small number NZ for tabulation to speed up the
        %    simulations
        tab = tabspace( layer, p );
        %  Green function table
        greentab = compgreentablayer( layer, tab );
        %  precompute Green function table
        %    for a more accurate simulation of the layer the number of
        %    wavelenghts should be increased
        greentab = parset( greentab, linspace( min(enei), max(enei), length(enei) ), op );
    end  
end
op.greentab = greentab;
p_p = gcp;
delete(p_p)
delete(gcp('nocreate'))
%% Green function  initialize BEM solver
bem = bemsolver( p, op );
%  initialize plane wave excitation
exc_flag = planewave( pol, dir, op );
%% Structure plot
f_tmp2 = figure;
%plot( p1, 'EdgeColor', 'b' );
plot(p, 'EdgeColor', 'r', 'FaceAlpha', 0.5)
%pos1 = quad(p);
%plot3( pos1( :, 1 ), pos1( :, 2 ), pos1( :, 3 ), 'r.' ,'MarkerSize',1 );
axis on
pbaspect([1 1 1])
title('Structure geometry')
%material dull 
%shading interp
%lighting none
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
saveas(f_tmp2,'data/str.fig')
%min( p.pos )
%% Field plot\
%%  computation of electric field
if field_set==1
    sig = bem \ exc_flag( p, enei );
    %  object for electric field
    %    MINDIST controls the minimal distance of the field points to the
    %    particle boundary, MESHFIELD must receive the OP structure which also
    %    stores the table of precomputed reflected Green functions
    emesh = meshfield( p, x, 0, z, op, 'mindist', 0.15, 'nmax', 2000 );
    %  induced and incoming electric field
    e = emesh( sig ) + emesh( exc_flag.field( emesh.pt, enei ) );
    %  norm of electric field
    ee = sqrt( dot( e, e, 3 ) );

    %%  plot electric field
    figure
    imagesc( x( : ), z( : ), log10(ee) ); 
    hold on;
    plot( [ min( x( : ) ), max( x( : ) ) ], [ 0, 0 ], 'w--' );
    plot( [ min( x( : ) ), max( x( : ) ) ], [ sub_t2, sub_t2 ], 'w--' );
    colorbar;  %colormap hot( 255 );
    axis([-inf inf -inf inf])
    daspect([1 1 1])
    xlabel( 'x (nm)' );
    ylabel( 'z (nm)' );
    title( 'Electric field XZ' );
    set( gca, 'YDir', 'norm' );
end
end




