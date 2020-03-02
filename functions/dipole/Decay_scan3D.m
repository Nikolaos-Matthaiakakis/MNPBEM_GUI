function Decay_scan3D(  d_length, exc, bem, ene2, enei2, p, imp)
    delete('data/LDOS_3D\*') % Delete old files
    multiWaitbar( 'LDOS 3D', 0, 'Color', 'g', 'CanCancel', 'on' );
    for ien=1:length(enei2)
        cancelRun() % Cancel function
        sig = bem \ exc(p, enei2( ien ) ); %  Surface charges 
        [ tot( ien, :, : ), rad( ien, :, : ) ] = exc.decayrate( sig ); %  Total and radiative decay rate       
        if length (ene2)<10
            figure
            subplot(2,1,1)
            semilogy( imp, sum(tot( ien, :, : ),3), '-') 
            xlabel( 'x (nm)' );
            ylabel( 'Total decay rate' );
            title(sprintf([ 'Total decay ' num2str(ene2(ien)) ' eV' ]));
            subplot(2,1,2)
            semilogy( imp, sum(rad( ien, :, : ),3), '-')
            xlabel( 'x (nm)' );
            ylabel( 'Radiative decay rate' );
            title(sprintf([ 'Radiative decay ' num2str(ene2(ien)) ' eV' ]));
            saveas(gcf,sprintf(['data/LDOS_3D/' sprintf( '%04d', ien ) 'LDOS3D_' num2str(ene2(ien)) 'eV' '.fig']))
            saveas(gcf,sprintf(['data/LDOS_3D/' sprintf( '%04d', ien ) 'LDOS3D_' num2str(ene2(ien)) 'eV' '.png']))      
            close             
        end
        multiWaitbar( 'LDOS 3D', ien / numel( enei2 ) );
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    if length (ene2)<10
        imagefiles = dir('data/LDOS_3D/*.png') ; % Find files in folder       
        nfiles = length(imagefiles); % Number of files found
        for ii=1:nfiles
            currentfilename = imagefiles(ii).name;
            currentimage = imread(sprintf(['data/LDOS_3D/' currentfilename]));
            images{ii} = currentimage; % Load image files
        end
        f=figure('Position', [100, 100, 1300, 800]);      
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
            name='LDOS 3D ';
            set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
            plotterfcn(vars_b,vars_im,ene2, name)  
        else
            imshow(images{1})
        end 
    else
        %  density plot of total scattering rate (LDOS)
        figure
        subplot(1,2,1)
        imagesc( ene2, imp, real(log10( sum( tot, 3 ) ) )  .' ); 
        hold on
        %  plot disk edge
        plot( ene2, 0 * ene2 + d_length/2, 'w--' );
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'LDOS 3D (Tot)' );
        subplot(1,2,2)
        imagesc( ene2, imp, real(log10( sum( rad, 3 ) ) ) .' ); 
        hold on
        %  plot disk edge
        plot( ene2, 0 * ene2 + d_length/2, 'w--' );
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'LDOS 3D (Rad)' );
        saveas(gcf,'data/LDOS_3D.fig')
        saveas(gcf,'data/LDOS_3D.png')
    end
end
