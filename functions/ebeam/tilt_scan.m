function tilt_scan(op, ene2, enei2, p, imp_spec, width, vel, spec_s)
    %[ psurf, pbulk ] = deal( zeros( numel( imp ), numel( enei ) )); %  Surface/bulk loss matrix init
    %CL_sca=deal( zeros( numel( 1:20:181 ), numel( enei ), 3 )); %  CL spectra matrix init
    multiWaitbar( 'Tilt scan', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths 
    misc.atomicunits;
    j=1;
    angle_r=360;
    step=10;
    for ang=1:step:angle_r+1
        p2 = rot( p, ang-1 , [1,0,0] ); %  Rotation around x (degrees) 
        bem1 = bemsolver( p2, op );  
        for ien = 1 : length( enei2 )            
            drawnow() % Reads new inputs
            cancelRun() % Cancel function              
            exc_CL = electronbeam( p2, imp_spec, width, vel, op );                                           
            sig = bem1 \ exc_CL( enei2( ien ) ); %  Surface charges        
            [ psurf( :, ien ), pbulk( :, ien ) ]= exc_CL.loss( sig ); %  EELS losses           
            ffang=farfield( spec_s,sig); % Farfield 
            [sca_ang,~] = scattering(ffang); %  Farfield scattering
            CL_sca(j,ien,1)=sca_ang.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei2(ien)); %  Farfield scattering spectra#
            %-------------------------% p pol
            ffp=ffang;
            ffp.e(:,2)=0; % E field in x direction set to 0
            ffp.h(:,1)=0; % B field in y direction set to 0
            [sca,~]=scattering(ffp); %  Farfield calculation
            CL_sca(j,ien,2)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei2(ien)); %  Farfield scattering spectra
            %-------------------------% s pol
            ffs=ffang;
            ffs.e(:,1)=0; % E field in y direction set to 0
            ffs.h(:,2)=0; % B field in x direction set to 0
            [sca,~]=scattering(ffs); %  Farfield calculation
            CL_sca(j,ien,3)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei2(ien)); %  Farfield scattering spectra
        end
        j=j+1;
        multiWaitbar( 'Tilt scan', j / numel( (1:step:angle_r+1)-1 ) );      
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    f=figure('Position', [100, 100, 600, 800]); 
    subplot(3,1,1);    
    imagesc( ene2, (1:step:angle_r+1)-1, CL_sca(:,:,1) );   
    hold on
    plot( ene2, 0 * ene2 + angle_r/2, 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Tilt scan' );
    hold off 
    
    subplot(3,1,2);
    imagesc( ene2, (1:step:angle_r+1)-1, CL_sca(:,:,2) );   
    hold on
    plot( ene2, 0 * ene2 + angle_r/2, 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Tilt scan p' );
    
    subplot(3,1,3);
    imagesc( ene2, (1:step:angle_r+1)-1, CL_sca(:,:,3) );   
    hold on
    plot( ene2, 0 * ene2 + angle_r/2, 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Tilt scan s' );
    
    saveas(gcf,'data/tilt_scan.fig')
    saveas(gcf,'data/tilt_scan.png')
end