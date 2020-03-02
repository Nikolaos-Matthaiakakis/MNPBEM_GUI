%%  L maps calculation, runs for impact instead of imp 
function L_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, p )
    name=["L", "L p", "L s", "L_a", "L_a p", "L_a s"];
    [ CL_sca_m_CL ] = deal( zeros( size( impact, 1 ), length( enei2 ) ) ); % Preloads matrix size
    CL_sca_m_CLp=CL_sca_m_CL; 
    CL_sca_m_CLs=CL_sca_m_CL;
    CL_sca_m_CLa=CL_sca_m_CL;
    CL_sca_m_CLap=CL_sca_m_CL;
    CL_sca_m_CLas=CL_sca_m_CL;
    multiWaitbar( 'Luminescence maps', 0, 'Color', 'g', 'CanCancel', 'on' );    
    for ien = 1 : length( enei2 ) %  Loop over energies
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL_M(p, enei2( ien ) ); %  Surface charges 
        % Whole sphere
        ff_m=farfield(spec_s,sig); % Farfield for whole sphere
        [sca_m,~] = scattering(ff_m); %  Farfield scattering for whole sphere
        CL_sca_m_CL(:,ien)=sca_m;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        %---- P pol
        ff_mp=ff_m;
        ff_mp.e(:,2,:)=0;%ff_m.e(:,1)=0; % E field in y direction set to 0
        ff_mp.h(:,1,:)=0; % B field in x direction set to 0
        [sca_mp,~] = scattering(ff_mp); %  Farfield scattering for whole sphere
        CL_sca_m_CLp(:,ien)=sca_mp;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        %---- S pol
        ff_ms=ff_m;
        ff_ms.e(:,1,:)=0; % E field in x direction set to 0
        ff_ms.h(:,2,:)=0; % B field in y direction set to 0
        [sca_ms,~] = scattering(ff_ms); %  Farfield scattering for whole sphere
        CL_sca_m_CLs(:,ien)=sca_ms;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
        % Angular
        ff_ma=farfield(spec,sig); % Farfield for angular  
        [sca_ma,~] = scattering(ff_ma); %  Farfield scattering for Angular
        CL_sca_m_CLa(:,ien)=sca_ma;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
        %---- P pol
        ff_map=ff_ma;
        ff_map.e(:,2,:)=0;%ff_m.e(:,1)=0; % E field in y direction set to 0
        ff_map.h(:,1,:)=0; % B field in x direction set to 0
        [sca_map,~] = scattering(ff_map); %  Farfield scattering for Angular
        CL_sca_m_CLap(:,ien)=sca_map;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
        %---- S pol
        ff_mas=ff_ma;
        ff_mas.e(:,1,:)=0; % E field in x direction set to 0
        ff_mas.h(:,2,:)=0; % B field in y direction set to 0
        [sca_mas,~] = scattering(ff_mas); %  Farfield scattering for Angular
        CL_sca_m_CLas(:,ien)=sca_mas;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
        multiWaitbar( 'Luminescence maps', ien / numel( enei2 ) );   
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar 
%%  Plot-section 2: L maps    
    figure; 
    delete('data/luminescence_map\*') % Delete old files    
    for ien = 1 : length(enei2)               
        prob = reshape( CL_sca_m_CL( :, ien ), size( xe ) ); %  Reshape CL spectra
        probp = reshape( CL_sca_m_CLp( :, ien ), size( xe ) );
        probs = reshape( CL_sca_m_CLs( :, ien ), size( xe ) );
        proba = reshape( CL_sca_m_CLa( :, ien ), size( xe ) );
        probap = reshape( CL_sca_m_CLap( :, ien ), size( xe ) );
        probas = reshape( CL_sca_m_CLas( :, ien ), size( xe ) );
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
        imagesc( xx, yy,  probs ); %  Add second part of structure
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
        saveas(gcf,['data/luminescence_map/' sprintf( '%04d', ien ) 'L_' num2str(ene2(ien)) 'eV' '.png'])
        saveas(gcf,['data/luminescence_map/' sprintf( '%04d', ien ) 'L_' num2str(ene2(ien)) 'eV' '.fig'])
    end 
    close
    imagefiles = dir('data/luminescence_map/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/luminescence_map/' currentfilename]));
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
        name='Luminescence ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end