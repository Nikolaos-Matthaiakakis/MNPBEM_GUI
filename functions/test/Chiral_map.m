%%  CL maps calculation, runs for impact instead of imp 
function Chiral_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en)
    name=["CL", "CL lcp", "CL rcp", "CL_a", "CL_a lcp", "CL_a rcp"];
    [ CL_sca_m_CL ] = deal( zeros( size( impact, 1 ), length( enei2 ) ) ); % Preloads matrix size
    CL_sca_m_CLlcp=CL_sca_m_CL; 
    CL_sca_m_CLrcp=CL_sca_m_CL;
    CL_sca_m_CLa=CL_sca_m_CL;
    CL_sca_m_CLalcp=CL_sca_m_CL;
    CL_sca_m_CLarcp=CL_sca_m_CL;
    epsilon_0=8.85*10^-14; %F/cm permittivity of vacuum
    c = 2.99792458*10^8; %m/s speed of light
    multiWaitbar( 'Chiral CL maps', 0, 'Color', 'g', 'CanCancel', 'on' );    
    for ien = 1 : length( enei2 ) %  Loop over energies
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL_M( enei2( ien ) ); %  Surface charges
        % Whole sphere
        ff_m=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_m,~] = scattering(ff_m); %  Farfield scattering for whole sphere
        CL_sca_m_CL(:,ien)=sca_m/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        %---- LCP pol   
        test=size(ff_m.e)
        ff_mlcp=ff_m;
        ff_mlcp.e(:,1,:)=abs(real(ff_m.e(:,1,:))-1i*imag(ff_m.e(:,2,:)));
        ff_mlcp.e(:,2,:)=0;
        ff_mlcp.h(:,2,:)=abs(real(ff_m.e(:,2,:))-1i*imag(ff_m.e(:,1,:)));
        ff_mlcp.h(:,2,:)=1;
        [sca_m,~] = scattering(ff_mlcp); %  Farfield scattering for whole sphere
        CL_sca_m_CLlcp(:,ien)=sca_m/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        %---- RCP pol
        ff_mrcp=ff_m;
        ff_mrcp.e(:,1,:)=abs(real(ff_m.e(:,1,:))+1i*imag(ff_m.e(:,2,:)));
        ff_mrcp.e(:,2,:)=0;
        ff_mrcp.h(:,2,:)=abs(real(ff_m.e(:,2,:))+1i*imag(ff_m.e(:,1,:)));
        ff_mrcp.h(:,2,:)=1;
        %ff_mrcp.h=epsilon_0*c*(abs(real(ff_m.h)+1i*imag(ff_m.h)).^2)/4;
        [sca_m,~] = scattering(ff_mrcp); %  Farfield scattering for whole sphere
        CL_sca_m_CLrcp(:,ien)=sca_m/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        % Angular
        ff_ma=farfield(spec,sig); % Farfield for angular  
        [sca_ma,~] = scattering(ff_ma); %  Farfield scattering for Angular
        CL_sca_m_CLa(:,ien)=sca_ma/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
        %---- LCP pol   
        ff_malcp=ff_ma;
        ff_malcp.e=epsilon_0*c*(abs(real(ff_ma.e)-1i*imag(ff_ma.e)).^2)/4;
        [sca_m,~] = scattering(ff_malcp); %  Farfield scattering for whole sphere
        CL_sca_m_CLalcp(:,ien)=sca_m/(2*pi/enei2(ien) ); %  Farfield scattering spectra 
        %---- RCP pol
        ff_marcp=ff_ma;
        ff_marcp.e=epsilon_0*c*(abs(real(ff_ma.e)+1i*imag(ff_ma.e)).^2)/4;
        [sca_m,~] = scattering(ff_marcp); %  Farfield scattering for whole sphere
        CL_sca_m_CLarcp(:,ien)=sca_m/(2*pi/enei2(ien) ); %  Farfield scattering spectra 
        multiWaitbar( 'Chiral CL maps', ien / numel( enei2 ) );   
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar 
%%  Plot-section 2: CL maps    
    figure; 
    delete('data/chiral_map\*') % Delete old files    
    for ien = 1 : length(enei2)               
        prob = reshape( CL_sca_m_CL( :, ien ), size( xe ) ); %  Reshape CL spectra
        probp = reshape( CL_sca_m_CLlcp( :, ien ), size( xe ) );
        probs = reshape( CL_sca_m_CLrcp( :, ien ), size( xe ) );
        proba = reshape( CL_sca_m_CLa( :, ien ), size( xe ) );
        probap = reshape( CL_sca_m_CLalcp( :, ien ), size( xe ) );
        probas = reshape( CL_sca_m_CLarcp( :, ien ), size( xe ) );
        subplot(2,3,1); %  Defines subplot window
        imagesc( xx, yy, prob ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(1));
        subplot(2,3,2); %  Defines subplot window
        imagesc( xx, yy, probp ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(2));
        subplot(2,3,3); %  Defines subplot window
        imagesc( xx, yy, probs ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(3));
        subplot(2,3,4); %  Defines subplot window
        imagesc( xx, yy, proba ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(4));
        subplot(2,3,5); %  Defines subplot window
        imagesc( xx, yy, probap ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(5));
        subplot(2,3,6); %  Defines subplot window
        imagesc( xx, yy, probas ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        title(name(6));
        %set( gca, 'YDir', 'norm' );
        saveas(gcf,['data/chiral_map/' sprintf( '%04d', ien ) 'chiral_' num2str(ene2(ien)) 'eV' '.png'])
        saveas(gcf,['data/chiral_map/' sprintf( '%04d', ien ) 'chiral_' num2str(ene2(ien)) 'eV' '.fig'])
    end 
    close
    imagefiles = dir('data/chiral_map/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/chiral_map/' currentfilename]));
        images{ii} = currentimage; % Load image files
    end
    f=figure('Position', [100, 100, 1000, 500]);     
    % Assign slider properties
    if length(ene2)>1
        val=1;
        b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'Value',val, 'Min',1, 'Max',length(ene2), 'SliderStep', [1/(length(ene2)-1) 1/(length(ene2)-1)]);
        bgcolor = f.Color;
        bl1 = uicontrol('Parent',f,'Style','text','Position',[50,24,23,23],...
                'String',min(ene2),'BackgroundColor',bgcolor);
        bl2 = uicontrol('Parent',f,'Style','text','Position',[500,24,23,23],...
                'String',max(ene2),'BackgroundColor',bgcolor);
        bl3 = uicontrol('Parent',f,'Style','text','Position',[240,5,100,13],...
                'String','eV','BackgroundColor',bgcolor);    
        % Set up callbacks
        vars_im=struct('images',images);
        vars_b=struct('b',b);
        name='Chiral ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end   