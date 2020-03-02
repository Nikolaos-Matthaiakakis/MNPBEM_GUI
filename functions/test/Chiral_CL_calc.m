%%  EELS loss and CL spectra calculation   
function Chiral_CL_calc(op, vel, exc, bem, Double_en, d_length, epstab, imp, imp_spec,ene, enei, p, spec, spec_s, exc_CL, typeStr)
    [ psurf, pbulk ] = deal( zeros( numel( imp ), numel( enei ) )); %  Surface/bulk loss matrix init
    CL_sca=deal( zeros( numel( imp ), numel( enei ), 3 )); %  CL spectra matrix init
    CL_sca_s=CL_sca;
    CL_sca_CL=CL_sca;
    CL_sca_s_CL=CL_sca;
    epsilon_0=8.85*10^-14; %F/cm permittivity of vacuum
    c = 2.99792458*10^8; %m/s speed of light
    multiWaitbar( 'EELS Loss/ CL', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths    
    for ien = 1 : length( enei )
        drawnow() % Reads new inputs
        cancelRun() % Cancel function        
        if length(imp)==1
            sig = bem \ exc_CL( enei( ien ) ); %  Surface charges 
            [ psurf( :, ien ), pbulk( :, ien ) ]= exc_CL.loss( sig ); %  EELS losses
        else
            sig = bem \ exc( enei( ien ) ); %  Surface charges 
            [ psurf( :, ien ), pbulk( :, ien ) ]= exc.loss( sig ); %  EELS losses
        end
        ff=farfield(spec,sig); % Farfield              
        %-------------------------% Both polarizations Angular
        [sca,~] = scattering(ff); %  Farfield scattering
        CL_sca(:,ien,1)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% p pol
        ffp=ff;
        ffp.e=epsilon_0*c*(abs(real(ff.e)-1i*imag(ff.e)).^2)/4;
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca(:,ien,2)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% s pol
        ffp=ff;
        ffp.e=epsilon_0*c*(abs(real(ff.e)+1i*imag(ff.e)).^2)/4;
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca(:,ien,3)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra       
        %-------------------------% Both polarizations whole sphere
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
        CL_sca_s(:,ien,1)=sca_s/(2*pi/enei(ien) ); %  Farfield scattering spectra for whole sphere
        
        
        %-------------------------% p pol
        ffp=ff_s;
        ffp.e=epsilon_0*c*(abs(real(ff_s.e)-1i*imag(ff_s.e)).^2)/4;
        [sca_s,~]=scattering(ffp); %  Farfield calculation
        CL_sca_s(:,ien,2)=sca_s/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% s pol
        ffp=ff_s;
        ffp.e=epsilon_0*c*(abs(real(ff_s.e)+1i*imag(ff_s.e)).^2)/4;
        [sca_s,~]=scattering(ffp); %  Farfield calculation
        CL_sca_s(:,ien,3)=sca_s/(2*pi/enei(ien) ); %  Farfield scattering spectra 
        %-------------------------%
        
        
        
        
        %-------------------------% repeats for imp_spec        
        sig = bem \ exc_CL( enei( ien ) ); %  Surface charges  
        ff=farfield(spec,sig); % Farfield
        %-------------------------% Both polarizations Angular
        [sca,~] = scattering(ff); %  Farfield scattering
        CL_sca_CL(:,ien,1)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% p pol
        
        
        cp=1;
        ffp=ff;
        ffs=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        %ffp.e=epsilon_0*c*(abs(real(ff.e(:,1,:))-1i*imag(ff.e(:,2,:))).^2)/4;
        [sca,~]=scattering_chiral(ffp,ffs,cp); %  Farfield calculation
        CL_sca_CL(:,ien,2)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra        
        
        
        %-------------------------% 
        cp=2;
        ffp=ff;
        ffs=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        ffp.e=epsilon_0*c*(abs(real(ff.e(:,1,:))-1i*imag(ff.e(:,2,:))).^2)/4;
        [sca,~]=scattering_chiral(ffp,ffs,cp); %  Farfield calculation
        CL_sca_CL(:,ien,3)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra 

        %-------------------------% Both polarizations whole sphere
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
        CL_sca_s_CL(:,ien,1)=sca_s/(2*pi/enei(ien) ); %  Farfield scattering spectra for whole sphere
        %-------------------------% p pol
        cp=1;
        ffp=ff_s;
        ffs=ff_s;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        %ffp.e=epsilon_0*c*(abs(real(ff.e(:,1,:))-1i*imag(ff.e(:,2,:))).^2)/4;
        [sca,~]=scattering_chiral(ffp,ffs,cp); %  Farfield calculation
        CL_sca_s_CL(:,ien,2)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra        
        
        
        %-------------------------% 
        cp=2;
        ffp=ff_s;
        ffs=ff_s;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        %ffp.e=epsilon_0*c*(abs(real(ff.e(:,1,:))-1i*imag(ff.e(:,2,:))).^2)/4;
        [sca,~]=scattering_chiral(ffp,ffs,cp); %  Farfield calculation
        CL_sca_s_CL(:,ien,3)=sca/(2*pi/enei(ien) ); %  Farfield scattering spectra

        multiWaitbar( 'EELS Loss/ CL', ien / numel( enei ) );
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar    

    %% Plot-section 1: Structure-mesh, EELS probability figure, CL spectra  
    figure; %  Plots EELS loss    
    %-------------------------% Plots CL spectra for specific spectral energy (whole sphere)
    subplot(2,2,1); %  Defines subplot window
    if length (imp)==1 %  Solution for only one beam position
        plot(p,'EdgeColor','b','nvec',1) % Run to verify face orientation
        if Double_en==0
            pbaspect([1 1 1])
            axis([-d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5])
        else
            pbaspect([2 1 1])
            axis_tmp=d_length+d_length/1.5+d_length/2.5+d_length/1.5+d_length/5;
            axis([-d_length-d_length/1.5-d_length/2.5 d_length/1.5+d_length/5 -axis_tmp/4 axis_tmp/4 -axis_tmp/4 axis_tmp/4])
        end
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        zlabel( 'z (nm)' );
    else %  Solution for many positions
        imagesc( ene, imp, CL_sca_s(:,:,1) );   
        hold on
        plot( ene, 0 * ene + d_length/2, 'w--' );
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'CL spectra' );
        hold off
    end 
    subplot(2,2,3)  
    if length (imp)==1 %  Solution for only one beam position
        if typeStr==0 %  Solution for sphere     
            mie = miesolver( epstab{ 2 }, epstab{ 1 }, d_length, op, 'lmax', 40 ); %  Mie solver
            if (abs(imp_spec(1))-d_length/2)>0
                plot( ene, psurf, 'o-',ene, pbulk,'.-',ene, psurf+pbulk,'-', ene, mie.loss( imp_spec(1)-d_length/2, enei, vel ),'*-') 
                legend( 'surface', 'bulk','EELS loss','Mie' );
            else
                plot( ene, psurf, 'o-',ene, pbulk,'.-',ene, psurf+pbulk,'-') 
                legend( 'surface', 'bulk','EELS loss' );
            end
        else %  Solution for other nanostructures
            plot( ene, psurf, 'o-',ene, pbulk,'.-',ene, psurf+pbulk,'-');
            legend( 'surface', 'bulk','EELS loss' );
        end
        %pbaspect([1 1 1])
        xlabel( 'Loss energy (eV)' );
        ylabel( 'Loss probability (eV^{-1})' );   
    else %  Solution for many positions
        imagesc( ene, imp, psurf + pbulk );  
        hold on
        plot( ene, 0 * ene + d_length/2, 'w--' ); %   Plots structure boundary 
        %pbaspect([1 1 1])
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Loss energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'Loss probability (eV^{-1})' );
        hold off
    end
    hold off  
    %-------------------------% Plots CL spectra for specific beam position (whole)
    subplot(2,2,2);
    plot(ene, CL_sca_s(1,:,2),'.-',ene, CL_sca_s(1,:,3),'-');  
    legend(  'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'CL spectra' );
    %-------------------------% Plots CL spectra for specific beam position (angular)
    subplot(2,2,4); %  Defines subplot window
    %if length (imp)==1 %  Solution for only one beam position
    plot( ene, CL_sca_CL(1,:,2),'.-',ene, CL_sca_CL(1,:,3),'-');
    legend( 'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'Angular CL spectra' );
    
    saveas(gcf,'data/spectra.fig')
    saveas(gcf,'data/spectra.png')
end