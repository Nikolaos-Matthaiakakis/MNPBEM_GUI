function Decay_map( exc, bem, ene2, enei2, p, Double_en, xe, xx, yy)
    delete('data/decay_maps\*') % Delete old files
    multiWaitbar( 'Total and radiative decay (map)', 0, 'Color', 'g', 'CanCancel', 'on' );
    for ien=1:length(enei2)
        cancelRun() % Cancel function
        sig = bem \ exc(p, enei2( ien ) ); %  Surface charges 
        [ tot(  :, ien ), rad(  :, ien ) ] = exc.decayrate( sig ); %  Total and radiative decay rate
    end
    %%  map plots
    for ien=1:length(enei2)   
        probt = reshape( tot( :, ien ), size( xe ) );
        probr = reshape( rad( :, ien ), size( xe ) );
        %  density plot of total scattering rate (LDOS)
        figure
        subplot(1,2,1)
        imagesc( xx, yy, real(probt) );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'x (nm)' );
        ylabel( 'y (nm)' );
        title( 'Total decay map' );
        subplot(1,2,2)
        imagesc( xx, yy, real(probr) );
        if Double_en==0
            pbaspect([1 1 1])
        else
            pbaspect([2 1 1])
        end
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'x (nm)' );
        ylabel( 'y (nm)' );
        title( 'Radiative decay map' );
        saveas(gcf,sprintf(['data/decay_maps/' sprintf( '%04d', ien ) 'Dec_map_' num2str(ene2(ien)) 'eV' '.fig']))
        saveas(gcf,sprintf(['data/decay_maps/' sprintf( '%04d', ien ) 'Dec_map_' num2str(ene2(ien)) 'eV' '.png']))      
        close             
        multiWaitbar( 'Total and radiative decay (map)', ien / numel( enei2 ) );
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    imagefiles = dir('data/decay_maps/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/decay_maps/' currentfilename]));
        images{ii} = currentimage; % Load image files
    end
    f=figure('Position', [100, 100, 1000, 800]);     
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
        name='Decay map ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end
