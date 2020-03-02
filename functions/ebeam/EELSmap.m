%%  EELS maps calculation, runs for impact instead of imp 
function EELSmap(impact, ene2, enei2, exc_EELS_M, xx, yy, bem, xe, Double_en)
    [ psurf, pbulk ] = deal( zeros( size( impact, 1 ), length( enei2 ) ) ); %  Reshape psurf and pbulk for EELS map
    multiWaitbar( 'EELS map', 0, 'Color', 'g', 'CanCancel', 'on' ); 
    for ien = 1 : length( enei2 ) %  Loop over energies
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_EELS_M( enei2( ien ) ); %  Surface charges        
        [ psurf( :, ien ), pbulk( :, ien ) ] = exc_EELS_M.loss( sig ); %  EELS losses
        multiWaitbar( 'EELS maps', ien / numel( enei2 ) );
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar
%%  Plot-section 4: EELS maps     
    figure; %  Plot EELS maps
    delete('data/eels_map\*') % Delete old files
    name='EELS ';
    for ien = 1 : length(enei2)
        prob = reshape( psurf( :, ien ) + pbulk( :, ien ), size( xe ) ); %  Reshape loss probability
        imagesc( xx, yy, prob ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        %colormap hot( 255 ); 
        colorbar;
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' )
        saveas(gcf,['data/eels_map/' sprintf( '%04d', ien ) 'EELS_' num2str(ene2(ien)) 'eV' '.png'])
        saveas(gcf,['data/eels_map/' sprintf( '%04d', ien ) 'EELS_' num2str(ene2(ien)) 'eV' '.fig'])         
    end    
    close
    imagefiles = dir('data/eels_map/*.png') ; % Find files        
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/eels_map/' currentfilename]));
        images{ii} = currentimage; % Load Image files
    end
    f=figure;
    if length(ene2)>1
        val=1;
        % Assign slider properties
        be = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'Value',val, 'Min',1, 'Max',length(ene2), 'SliderStep', [1/(length(ene2)-1) 1/(length(ene2)-1)]);
        bgecolor = f.Color;
        bel1 = uicontrol('Parent',f,'Style','text','Position',[50,24,23,23],...
                'String',min(ene2),'BackgroundColor',bgecolor);
        bel2 = uicontrol('Parent',f,'Style','text','Position',[500,24,23,23],...
                'String',max(ene2),'BackgroundColor',bgecolor);
        bel3 = uicontrol('Parent',f,'Style','text','Position',[240,5,100,13],...
                'String','eV','BackgroundColor',bgecolor);    
        % Set up callbacks
        varse_im=struct('images',images);
        varse_b=struct('b',be);
        set(be,'Callback',{@slider1_callback,varse_b,varse_im,ene2, name});
        plotterfcn(varse_b,varse_im,ene2, name)  
    else
        imshow(images{1})
    end
end