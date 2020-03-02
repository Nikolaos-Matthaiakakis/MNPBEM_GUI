%%  EELS loss and CL spectra calculation   
function FF_sca_calc_den(exc_CL, bem, Double_en, d_length, ene, enei, p, spec, spec_s, shift_l)
    multiWaitbar( 'Decay rate and far-field spectra', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths 
    CL_sca=deal( zeros( 1, numel( enei ), 3 )); %  CL spectra matrix init
    CL_sca_s=CL_sca;
    scatter=deal( zeros( 1, numel( enei ) ));
    for j=1:3
        for ien = 1 : length( enei )
            drawnow() % Reads new inputs
            cancelRun() % Cancel function
            sig = bem \ exc_CL(p, enei( ien ) ); %  Surface charges        
            scatter(:, ien) = exc_CL.decayrate( sig );
            if j==1
                ff=farfield( exc_CL, spec, enei(ien) ); % Farfield    
            elseif j==2
                ff=farfield(spec,sig); % Farfield 
            else
                ff=farfield(spec,sig)+ farfield( exc_CL, spec, enei( ien ) ); % Farfield 
            end
            %-------------------------% Both polarizations Angular
            [sca,~] = scattering(ff); %  Farfield scattering
            CL_sca(ien,1,j)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
            %-------------------------% p pol
            ffp=ff;
            ffp.e(:,2)=0; % E field in x direction set to 0
            ffp.h(:,1)=0; % B field in y direction set to 0
            [sca,~]=scattering(ffp); %  Farfield calculation
            CL_sca(ien,2,j)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
            %-------------------------% s pol
            ffs=ff;
            ffs.e(:,1)=0; % E field in y direction set to 0
            ffs.h(:,2)=0; % B field in x direction set to 0
            [sca,~]=scattering(ffs); %  Farfield calculation
            CL_sca(ien,3,j)=sca;%/(2*pi/enei(ien) ); %  Farfield scattering spectra       
            %-------------------------% Both polarizations whole sphere
            if j==1
                ff_s=farfield( exc_CL, spec_s, enei(ien)  ); % Farfield    
            elseif j==2
                ff_s=farfield(spec_s,sig); % Farfield 
            else
                ff_s=farfield(spec_s,sig)+farfield( exc_CL, spec_s, enei( ien ) ); % Farfield 
            end
            [sca_s,~] = scattering(ff_s); %  Farfield scattering for whole sphere
            CL_sca_s(ien,1,j)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra for whole sphere
            %-------------------------% p pol
            ffp_s=ff_s;
            ffp_s.e(:,2)=0; % E field in x direction set to 0
            ffp_s.h(:,1)=0; % B field in y direction set to 0
            [sca_s,~]=scattering(ffp_s); %  Farfield calculation
            CL_sca_s(ien,2,j)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra
            %-------------------------% s pol
            ffs_s=ff_s;
            ffs_s.e(:,1)=0; % E field in y direction set to 0
            ffs_s.h(:,2)=0; % B field in x direction set to 0
            [sca_s,~]=scattering(ffs_s); %  Farfield calculation
            CL_sca_s(ien,3,j)=sca_s;%/(2*pi/enei(ien) ); %  Farfield scattering spectra 
            multiWaitbar( 'Decay rate and far-field spectra', ien / numel( enei ));
            %-------------------------% repeats for imp_spec        
        end
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
    plot( ene, CL_sca_s(:,1,1), 'o-',ene, CL_sca_s(:,1,2),'.-',ene, CL_sca_s(:,1,3),'-');  
    %legend( 'Decay rate');
    %pbaspect([1 1 1])
    legend( 'Dipole', 'Structure','Dip+Struct');
    xlabel( 'Energy (eV)' );
    ylabel( 'Scattering (a.u)' );   
    title( 'Far-field spectra' );
    hold off  
    %-------------------------% Plots CL spectra for specific beam position (whole)
    subplot(2,2,2);
    plot( ene, CL_sca_s(:,1,3)./CL_sca_s(:,1,1));
    legend( 'Total');
    xlabel( 'Energy (eV)' );
    ylabel( 'L/Ld' );
    title( 'Far-field spectra enhancement' );
    %-------------------------% Plots CL spectra for specific beam position (angular)
    subplot(2,2,4); %  Defines subplot window
    %if length (imp)==1 %  Solution for only one beam position
    plot( ene, CL_sca(:,1,3)./CL_sca(:,1,1));  
    legend( 'Total');
    xlabel( 'Energy (eV)' );
    ylabel( 'L/Ld' );
    title( 'Directional far-field spectra enhancement' );
    saveas(gcf,'data/spectra.fig')
    saveas(gcf,'data/spectra.png')
end