%%  EELS loss and CL spectra calculation   
function EELS_CL_calc(op, vel, exc, bem, Double_en, d_length, epstab, imp, imp_spec,ene, enei, p, spec, spec_s, exc_CL, typeStr, shift_l)
    [ psurf, pbulk ] = deal( zeros( numel( imp ), numel( enei ) )); %  Surface/bulk loss matrix init
    CL_sca=deal( zeros( numel( imp ), numel( enei ), 3 )); %  CL spectra matrix init
    CL_sca_s=CL_sca;
    CL_sca_CL=CL_sca;
    CL_sca_s_CL=CL_sca;
    misc.atomicunits;
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
        CL_sca(:,ien,1)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% p pol
        ffp=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca(:,ien,2)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs=ff;
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        [sca,~]=scattering(ffs); %  Farfield calculation
        CL_sca(:,ien,3)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra       
        %-------------------------% Both polarizations whole sphere
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
        CL_sca_s(:,ien,1)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra for whole sphere
        %-------------------------% p pol
        ffp_s=ff_s;
        ffp_s.e(:,2)=0; % E field in x direction set to 0
        ffp_s.h(:,1)=0; % B field in y direction set to 0
        [sca_s,~]=scattering(ffp_s); %  Farfield calculation
        CL_sca_s(:,ien,2)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs_s=ff_s;
        ffs_s.e(:,1)=0; % E field in y direction set to 0
        ffs_s.h(:,2)=0; % B field in x direction set to 0
        [sca_s,~]=scattering(ffs_s); %  Farfield calculation
        CL_sca_s(:,ien,3)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra 
        %-------------------------%
        %-------------------------% repeats for imp_spec        
        sig = bem \ exc_CL( enei( ien ) ); %  Surface charges  
        ff=farfield(spec,sig); % Farfield
        %-------------------------% Both polarizations Angular
        [sca,~] = scattering(ff); %  Farfield scattering
        CL_sca_CL(:,ien,1)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% p pol
        ffp=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca_CL(:,ien,2)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs=ff;
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        [sca,~]=scattering(ffs); %  Farfield calculation
        CL_sca_CL(:,ien,3)=sca.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% Both polarizations whole sphere
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
        CL_sca_s_CL(:,ien,1)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra for whole sphere
        %-------------------------% p pol
        ffp_s=ff_s;
        ffp_s.e(:,2)=0; % E field in x direction set to 0
        ffp_s.h(:,1)=0; % B field in y direction set to 0
        [sca_s,~]=scattering(ffp_s); %  Farfield calculation
        CL_sca_s_CL(:,ien,2)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs_s=ff_s;
        ffs_s.e(:,1)=0; % E field in y direction set to 0
        ffs_s.h(:,2)=0; % B field in x direction set to 0
        [sca_s,~]=scattering(ffs_s); %  Farfield calculation
        CL_sca_s_CL(:,ien,3)=sca_s.*fine ^ 2 / ( bohr * hartree * pi * vel )./(2*(pi^2)./enei(ien)); %  Farfield scattering spectra 
        %-------------------------% ellipticity for whole sphere
        fv_s=[mean(ff_s.e(:,1));mean(ff_s.e(:,2))];
        [delta_circ_ma_s(ien),n_circ_ma_s(ien),ar_s(ien),rs_s(ien,:)]=polellip(fv_s);
        %-------------------------% ellipticity for angular
        fv=[mean(ff.e(:,1));mean(ff.e(:,2))];
        [delta_circ_ma(ien),n_circ_ma(ien),ar(ien),rs(ien,:)]=polellip(fv); 
        % Uncomment for phase
        %ExPhi = angle(mean(ff_s.e(:,1)));
        %EyPhi = angle(mean(ff_s.e(:,2)));
        %phi(ien) = EyPhi-ExPhi;       
        multiWaitbar( 'EELS Loss/ CL', ien / numel( enei ) );
    end   
    % Uncomment for phase
    %figure
    %plot(ene, rad2deg(phi))   
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
            axis_tmp=(2*d_length+shift_l)/1.6;
            axis([-axis_tmp axis_tmp -axis_tmp/2 axis_tmp/2 -axis_tmp/2 axis_tmp/2])
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
        title( 'CL probability (eV^{-1})' );
        hold off
    end 
    subplot(2,2,3)  
    if length (imp)==1 %  Solution for only one beam position
        if typeStr==0 %  Solution for sphere     
            mie = miesolver( epstab{ 2 }, epstab{ 1 }, d_length, op, 'lmax', 40 ); %  Mie solver
            if (abs(imp_spec(1))-d_length/2)>0
                plot( ene, psurf, 'o-',ene, pbulk,'.-',ene, psurf+pbulk,'-', ene, mie.loss( abs(imp_spec(1))-d_length/2, enei, vel ),'*-') 
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
    plot( ene, CL_sca_s_CL(1,:,1), 'o-',ene, CL_sca_s_CL(1,:,2),'.-',ene, CL_sca_s_CL(1,:,3),'-');  
    legend( 'Total', 'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'CL probability (eV^{-1})' );
    %-------------------------% Plots CL spectra for specific beam position (angular)
    subplot(2,2,4); %  Defines subplot window
    %if length (imp)==1 %  Solution for only one beam position
    plot( ene, CL_sca_CL(1,:,1), 'o-',ene, CL_sca_CL(1,:,2),'.-',ene, CL_sca_CL(1,:,3),'-');
    legend( 'Total', 'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'Directional CL probability (eV^{-1})' );
    
    saveas(gcf,'data/spectra.fig')
    saveas(gcf,'data/spectra.png')
    %-------------------------% Plots ellipticity spectra for specific beam position (whole)
    figure
    subplot(2,1,1); %  Defines subplot window
    if length(enei)==1
        polellip(fv_s)
        %axis equal
        title({sprintf('Polarization Ellipse %1.1f eV',ene);sprintf('n=%1.1f theta=%1.1f %s', n_circ_ma_s , delta_circ_ma_s, rs_s{1})})
    else
        plot( ene, n_circ_ma_s, 'o-',ene, delta_circ_ma_s,'.-');
        legend( 'n','theta');
        xlabel( 'Energy (eV)' );
        ylabel( 'n, theta (Degress)' );
        title( 'Ellipticity spectra' );
    end
    %-------------------------% Plots ellipticity spectra for specific beam position (angular)
    subplot(2,1,2); %  Defines subplot window
    if length(enei)==1
        polellip(fv)
        %axis equal
        title({sprintf('Angular Polarization Ellipse %1.1f eV',ene);sprintf('n=%1.1f theta=%1.1f %s', n_circ_ma , delta_circ_ma, rs{1})})
    else
        plot( ene, n_circ_ma, 'o-',ene, delta_circ_ma,'.-');
        legend( 'n','theta');
        xlabel( 'Energy (eV)' );
        ylabel( 'n, theta (Degress)' );
        title( 'Directional ellipticity spectra' );
    end
    saveas(gcf,'data/ellipticity.fig')
    saveas(gcf,'data/ellipticity.png')
end