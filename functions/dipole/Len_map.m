%%  L maps calculation, runs for impact instead of imp 
function Len_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, p)
    name=["L/Ld", "L/Ld p", "L/Ld s", "L_a/Ld_a", "L_a/Ld_a p", "L_a/Ld_a s"];
    [ CL_sca_m_CL ] = deal( zeros( size( impact, 1 ), length( enei2 ) ) ); % Preloads matrix size
    CL_sca_m_CLp=CL_sca_m_CL; 
    CL_sca_m_CLs=CL_sca_m_CL;
    CL_sca_m_CLa=CL_sca_m_CL;
    CL_sca_m_CLap=CL_sca_m_CL;
    CL_sca_m_CLas=CL_sca_m_CL;    
    for j=1:2
        multiWaitbar( 'Luminescence enhancement maps', 0, 'Color', 'g', 'CanCancel', 'on' ); 
        for ien = 1 : length( enei2 ) %  Loop over energies
            drawnow() % Reads new inputs
            cancelRun() % Cancel function
            sig = bem \ exc_CL_M(p, enei2( ien ) ); %  Surface charges 
            % Whole sphere
            if j==1
                    ff_m=farfield( exc_CL_M, spec_s, enei2(ien)  ); % Farfield    
            else
                    ff_m=farfield(spec_s,sig)+farfield( exc_CL_M, spec_s, enei2( ien ) ); % Farfield 
            end
            [sca_m,~] = scattering(ff_m); %  Farfield scattering for whole sphere
            CL_sca_m_CL(:,ien,j)=sca_m;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
            %---- P pol
            ff_mp=ff_m;
            ff_mp.e(:,2,:)=0;%ff_m.e(:,1)=0; % E field in y direction set to 0
            ff_mp.h(:,1,:)=0; % B field in x direction set to 0
            [sca_mp,~] = scattering(ff_mp); %  Farfield scattering for whole sphere
            CL_sca_m_CLp(:,ien,j)=sca_mp;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
            %---- S pol
            ff_ms=ff_m;
            ff_ms.e(:,1,:)=0; % E field in x direction set to 0
            ff_ms.h(:,2,:)=0; % B field in y direction set to 0
            [sca_ms,~] = scattering(ff_ms); %  Farfield scattering for whole sphere
            CL_sca_m_CLs(:,ien,j)=sca_ms;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for whole sphere
            % Angular
            if j==1
                    ff_ma=farfield( exc_CL_M, spec, enei2(ien) ); % Farfield    
            else
                    ff_ma=farfield(spec,sig)+ farfield( exc_CL_M, spec, enei2( ien ) ); % Farfield 
            end
            [sca_ma,~] = scattering(ff_ma); %  Farfield scattering for Angular
            CL_sca_m_CLa(:,ien,j)=sca_ma;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
            %---- P pol
            ff_map=ff_ma;
            ff_map.e(:,2,:)=0;%ff_m.e(:,1)=0; % E field in y direction set to 0
            ff_map.h(:,1,:)=0; % B field in x direction set to 0
            [sca_map,~] = scattering(ff_map); %  Farfield scattering for Angular
            CL_sca_m_CLap(:,ien,j)=sca_map;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
            %---- S pol
            ff_mas=ff_ma;
            ff_mas.e(:,1,:)=0; % E field in x direction set to 0
            ff_mas.h(:,2,:)=0; % B field in y direction set to 0
            [sca_mas,~] = scattering(ff_mas); %  Farfield scattering for Angular
            CL_sca_m_CLas(:,ien,j)=sca_mas;%/(2*pi/enei2(ien) ); %  Farfield scattering spectra for Angular
            multiWaitbar( 'Luminescence enhancement maps', ien / numel( enei2 ) );   
        end    
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar 
%%  Plot-section 2: L maps    
    figure('Position', [100, 100, 1000, 500]);  
    delete('data/lum_enhance_map\*') % Delete old files    
    for ien = 1 : length(enei2)   
        for j=1:2
            prob(:,:,j) = reshape( CL_sca_m_CL( :, ien,j ), size( xe ) ); %  Reshape CL spectra
            probp(:,:,j) = reshape( CL_sca_m_CLp( :, ien,j ), size( xe ) );
            probs(:,:,j) = reshape( CL_sca_m_CLs( :, ien,j ), size( xe ) );
            proba(:,:,j) = reshape( CL_sca_m_CLa( :, ien,j ), size( xe ) );
            probap(:,:,j) = reshape( CL_sca_m_CLap( :, ien,j ), size( xe ) );
            probas(:,:,j) = reshape( CL_sca_m_CLas( :, ien,j ), size( xe ) );
        end
        subplot(2,3,1); %  Defines subplot window
        imagesc( xx, yy, prob(:,:,2)./prob(:,:,1) ); %  Add second part of structure
        colorbar
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
        imagesc( xx, yy, probp(:,:,2)./probp(:,:,1) ); %  Add second part of structure
        colorbar
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
        imagesc( xx, yy,  probs(:,:,2)./probs(:,:,1) ); %  Add second part of structure
        colorbar
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
        imagesc( xx, yy, proba(:,:,2)./proba(:,:,1) ); %  Add second part of structure
        colorbar
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
        imagesc( xx, yy, probap(:,:,2)./probap(:,:,1) ); %  Add second part of structure
        colorbar
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
        imagesc( xx, yy, probas(:,:,2)./probas(:,:,1) ); %  Add second part of structure
        colorbar
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
        saveas(gcf,['data/lum_enhance_map/' sprintf( '%04d', ien ) 'L_' num2str(ene2(ien)) 'eV' '.png'])
        saveas(gcf,['data/lum_enhance_map/' sprintf( '%04d', ien ) 'L_' num2str(ene2(ien)) 'eV' '.fig'])
    end 
    close
    imagefiles = dir('data/lum_enhance_map/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/lum_enhance_map/' currentfilename]));
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
        name='Luminescence enhancement ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end