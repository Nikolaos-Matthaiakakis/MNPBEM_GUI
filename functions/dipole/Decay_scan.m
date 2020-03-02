function Decay_scan( op, epstab, d_length, typeStr, exc, bem, ene2, enei2, dip_p, dip_d, p, imp)
    delete('data/decay_scan\*') % Delete old files
    multiWaitbar( 'Total and radiative decay (scan)', 0, 'Color', 'g', 'CanCancel', 'on' );
    if typeStr==0 %  Solution for sphere
        mie = miesolver( epstab{ 2 }, epstab{ 1 }, d_length, op, 'lmax', 40 ); %  Mie solver
    end
    for ien=1:length(enei2)
        cancelRun() % Cancel function
        sig = bem \ exc(p, enei2( ien ) ); %  Surface charges 
        [ tot( ien, : ), rad( ien, : ) ] = exc.decayrate( sig ); %  Total and radiative decay rate
        if typeStr==0 %  Solution for sphere
            [ tot0, rad0] = mie.decayrate( enei2( ien ),imp );
            if dip_d==1
                tot0p( ien, : )=tot0( :, 2 );
                rad0p( ien, : )=rad0( :, 2 );
            elseif dip_d==3
                tot0p( ien, : )=tot0( :, 1 );
                rad0p( ien, : )=rad0( :, 1 );
            end
        end      
        %%  line plots
        if length (ene2)<10
            figure
            subplot(2,1,1)
            if dip_d~=2 && typeStr==0 %  Solution for sphere
                semilogy( imp, tot( ien, : ), '-', imp, tot0p( ien, : ), 'g*'); 
                legend('Total-simulation','Total-Mie')
            else
                semilogy( imp, tot( ien, : ), '-')
            end    
            xlabel( 'x (nm)' );
            ylabel( 'Total decay rate' );
            title(sprintf([ 'Total decay ' num2str(ene2(ien)) ' eV' ]));
            subplot(2,1,2)
             if dip_d~=2 && typeStr==0 %  Solution for sphere
                semilogy( imp, rad( ien, : ), '-', imp, rad0p( ien, : ), 'g*'); 
                legend('Radiative-simulation','Radiative-Mie')
            else
                semilogy( imp, rad( ien, : ), '-')
            end    
            xlabel( 'x (nm)' );
            ylabel( 'Radiative decay rate' );
            title(sprintf([ 'Radiative decay ' num2str(ene2(ien)) ' eV' ]));
            saveas(gcf,sprintf(['data/decay_scan/' sprintf( '%04d', ien ) 'Dec_scan_' num2str(ene2(ien)) 'eV' '.fig']))
            saveas(gcf,sprintf(['data/decay_scan/' sprintf( '%04d', ien ) 'Dec_scan_' num2str(ene2(ien)) 'eV' '.png']))      
            close             
        end
        multiWaitbar( 'Total and radiative decay (scan)', ien / numel( enei2 ) );
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar  
    if length (ene2)<10
        imagefiles = dir('data/decay_scan/*.png') ; % Find files in folder       
        nfiles = length(imagefiles); % Number of files found
        for ii=1:nfiles
            currentfilename = imagefiles(ii).name;
            currentimage = imread(sprintf(['data/decay_scan/' currentfilename]));
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
            name='Decay scan ';
            set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
            plotterfcn(vars_b,vars_im,ene2, name)  
        else
            imshow(images{1})
        end 
    else
        %  density plot of total scattering rate (LDOS)
        figure
        subplot(1,2,1)
        imagesc( ene2, imp, real(log10(tot))  .' ); 
        hold on
        %  plot disk edge
        plot( ene2, 0 * ene2 + d_length/2, 'w--' );
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'LDOS (Tot)' );
        subplot(1,2,2)
        imagesc( ene2, imp, real(log10( rad)) .' ); 
        hold on
        %  plot disk edge
        plot( ene2, 0 * ene2 + d_length/2, 'w--' );
        set( gca, 'YDir', 'norm' );
        colorbar
        xlabel( 'Energy (eV)' );
        ylabel( 'x (nm)' );
        title( 'LDOS (Rad)' );
        saveas(gcf,'data/LDOS.fig')
        saveas(gcf,'data/LDOS.png')
    end
end
