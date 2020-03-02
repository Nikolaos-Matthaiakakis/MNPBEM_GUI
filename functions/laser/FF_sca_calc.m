%%  EELS loss and CL spectra calculation   
function FF_sca_calc(op, exc_CL, bem, Double_en, d_length, epstab, ene, enei, p, spec, spec_s, typeStr, shift_l)
    multiWaitbar( 'Scattering and far-field spectra', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths  
    CL_sca=deal( zeros( 1, numel( enei ), 3 )); %  CL spectra matrix init
    CL_sca_s=CL_sca;
    scatter=deal( zeros( 1, numel( enei ) ));
    for ien = 1 : length( enei )
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL(p, enei( ien ) ); %  Surface charges        
        scatter(:, ien) = exc_CL.sca( sig );
        extinct( :, ien ) = exc_CL.ext( sig );
        absorb( :, ien ) = exc_CL.abs( sig );
        ff=farfield(spec,sig); % Farfield              
        %-------------------------% Both polarizations Angular
        [sca,~] = scattering(ff); %  Farfield scattering
        CL_sca(:,ien,1)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% p pol
        ffp=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        [sca,~]=scattering(ffp); %  Farfield calculation
        CL_sca(:,ien,2)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs=ff;
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        [sca,~]=scattering(ffs); %  Farfield calculation
        CL_sca(:,ien,3)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra       
        %-------------------------% Both polarizations whole sphere
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
        CL_sca_s(:,ien,1)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra for whole sphere
        %-------------------------% p pol
        ffp_s=ff_s;
        ffp_s.e(:,2)=0; % E field in x direction set to 0
        ffp_s.h(:,1)=0; % B field in y direction set to 0
        [sca_s,~]=scattering(ffp_s); %  Farfield calculation
        CL_sca_s(:,ien,2)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
        %-------------------------% s pol
        ffs_s=ff_s;
        ffs_s.e(:,1)=0; % E field in y direction set to 0
        ffs_s.h(:,2)=0; % B field in x direction set to 0
        [sca_s,~]=scattering(ffs_s); %  Farfield calculation
        CL_sca_s(:,ien,3)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra 
        %-------------------------% repeats for imp_spec 
        %-------------------------% ellipticity for whole sphere
        fv_s=[mean(ff_s.e(:,1));mean(ff_s.e(:,2))];
        [delta_circ_ma_s(ien),n_circ_ma_s(ien),ar_s(ien),rs_s(ien,:)]=polellip(fv_s);
        %-------------------------% ellipticity for angular
        fv=[mean(ff.e(:,1));mean(ff.e(:,2))];
        [delta_circ_ma(ien),n_circ_ma(ien),ar(ien),rs(ien,:)]=polellip(fv);
        multiWaitbar( 'Scattering and far-field spectra', ien / numel( enei ) ); 
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar    

    %% Plot-section 1: Structure-mesh, EELS probability figure, CL spectra  
    figure; %  Plots EELS loss    
    %-------------------------% Plots CL spectra for specific spectral energy (whole sphere)
    subplot(2,2,1); %  Defines subplot window
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
    subplot(2,2,3)
    if typeStr==0 %  Solution for sphere         
        mie = miesolver( epstab{ 2 }, epstab{ 1 }, d_length, op, 'lmax', 40 ); %  Mie solver
        %plot( ene, scatter, 'o-',ene, extinct,'.-', ene, mie.sca( enei ),'*-')
        plot( ene, scatter, 'o-', ene, mie.sca( enei ),'*-') 
        %plot( ene,absorb);
        legend( 'Scattering','Mie' );
    else %  Solution for other nanostructures
        plot( ene, scatter, 'o-');%,ene, extinct,'.-');
        %legend( 'Scattering', 'extinction' );
        %plot( ene,absorb);
        %plot( ene,extinct, ene,absorb,'.-');  
        %plot( ene,scatter);
        legend( 'Scattering');%,'Exctinction' );
    end
    %pbaspect([1 1 1])
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering cross section (nm^2)' ); 
    %ylabel( 'Absorption cross section' ); 
    hold off  
    %-------------------------% Plots CL spectra for specific beam position (whole)
    subplot(2,2,2);
    plot( ene, CL_sca_s(1,:,1), 'o-',ene, CL_sca_s(1,:,2),'.-',ene, CL_sca_s(1,:,3),'-');  
    legend( 'Total', 'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'Far-field spectra' );
    %-------------------------% Plots CL spectra for specific beam position (angular)
    subplot(2,2,4); %  Defines subplot window
    %if length (imp)==1 %  Solution for only one beam position
    plot( ene, CL_sca(1,:,1), 'o-',ene, CL_sca(1,:,2),'.-',ene, CL_sca(1,:,3),'-');
    legend( 'Total', 'p','s');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );
    title( 'Directional far-field spectra' );
    saveas(gcf,'data/spectra.fig')
    saveas(gcf,'data/spectra.png')
    %-------------------------% Absorption and exctinction
    figure
    plot( ene,extinct, ene,absorb,'.-');  
    legend( 'Extinction','Absorption' );
    xlabel( 'Energy (eV)' );
    ylabel( 'Cross section (nm^2)' ); 
    title('Extinction and absorption cross section')
    saveas(gcf,'data/spectra_ext_abs.fig')
    saveas(gcf,'data/spectra_ext_abs.png')
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