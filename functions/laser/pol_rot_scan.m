function pol_rot_scan(p, ene2, enei2, spec_s, exc, bem)
    multiWaitbar( 'Polarization rotation scan', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths 
    step=5;        
    for ien = 1 : length( enei2 )            
        drawnow() % Reads new inputs
        cancelRun() % Cancel function                    

        sig = bem \ exc(p, enei2( ien ) ); %  Surface charges               
        ffang=farfield( spec_s,sig); % Farfield 
        [sca_ang,~] = scattering(ffang); %  Farfield scattering
        CL_sca(:,ien,1)=sca_ang;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra#
        %-------------------------% p pol
        ffp=ffang;
        ffp.e(:,2,:)=0; % E field in x direction set to 0
        ffp.h(:,1,:)=0; % B field in y direction set to 0
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca(:,ien,2)=sca;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs=ffang;
        ffs.e(:,1,:)=0; % E field in y direction set to 0
        ffs.h(:,2,:)=0; % B field in x direction set to 0
        [sca,~]=scattering(ffs); %  Farfield calculation
        CL_sca(:,ien,3)=sca;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra
    multiWaitbar( 'Polarization rotation scan', ien / numel( enei2 ) );      
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    f=figure('Position', [100, 100, 600, 800]); 
    subplot(3,1,1);    
    imagesc( ene2, (-90:step:90), CL_sca(:,:,1) );   
    hold on
    plot( ene2, 0 * ene2 , 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Polarization rotation scan' );
    hold off 
    
    subplot(3,1,2);
    imagesc( ene2, (-90:step:90), CL_sca(:,:,2) );   
    hold on
    plot( ene2, 0 * ene2 , 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Polarization rotation p' );
    
    subplot(3,1,3);
    imagesc( ene2, (-90:step:90), CL_sca(:,:,3) );   
    hold on
    plot( ene2, 0 * ene2 , 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Polarization rotation s' );
    
    saveas(gcf,'data/pol_rot.fig')
    saveas(gcf,'data/pol_rot.png')
end