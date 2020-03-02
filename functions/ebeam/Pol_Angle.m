function Pol_Angle(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, d_angle1)
    multiWaitbar( 'Ellipticity angle', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths 
    j=1;
    angle_r=180;
    step=4;
    for ang=1:step:angle_r+1
        for ien = 1 : length( enei2 )            
            drawnow() % Reads new inputs
            cancelRun() % Cancel function
            sig = bem \ exc_CL( enei2( ien ) ); %  Surface charges        
            [ psurf( :, ien ), pbulk( :, ien ) ]= exc_CL.loss( sig ); %  EELS losses
            angC=degtorad(angle_r-(ang-1)); % XZ angle in rad
            [spec_ang,~] = detInit(op,BEM_op,d_angle1,d_range1,angC,detQ,detRad);
            ffang=farfield(spec_ang,sig); % Farfield                        
            fv=[mean(ffang.e(:,1));mean(ffang.e(:,2))];
            [delta_circ_ma(j,ien),n_circ_ma(j,ien)]=polellip(fv);
        end
        j=j+1;
        multiWaitbar( 'Ellipticity angle', j / numel( (1:step:angle_r+1)-1 ) );      
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    f=figure('Position', [100, 100, 1000, 500]); 
    subplot(1,2,1);    
    imagesc( ene2, (1:step:angle_r+1)-1,  n_circ_ma(:,:) );   
    hold on
    plot( ene2, 0 * ene2 + angle_r/2, 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'n' );
    hold off 
    
    subplot(1,2,2);
    imagesc( ene2, (1:step:angle_r+1)-1, delta_circ_ma(:,:) );   
    hold on
    plot( ene2, 0 * ene2 + angle_r/2, 'w--' );
    set( gca, 'YDir', 'norm' );
    colorbar
    xlabel( 'Energy (eV)' );
    ylabel( 'Angle (Degrees)' );
    title( 'Theta' );   
     
    saveas(gcf,'data/Elliptangle.fig')
    saveas(gcf,'data/Elliptangle.png')
end