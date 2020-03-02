%%  CL maps calculation, runs for impact instead of imp  
% Reference(Imaging Chirality of optical fieldsnear achiral metal nanosctructures excited with linearly polarized light)
function Pol_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en)
    name=["n", "theta", "Directional n", "Directional theta"];
    [ n_circ_m ] = deal( zeros( size( impact, 1 ), length( enei2 ) ) ); % Preloads matrix size
    delta_circ_m=n_circ_m; 
    n_circ_ma=n_circ_m; 
    delta_circ_ma=n_circ_m; 
    multiWaitbar( 'Ellipticity maps', 0, 'Color', 'g', 'CanCancel', 'on' );    
    for ien = 1 : length( enei2 ) %  Loop over energies
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL_M( enei2( ien ) ); %  Surface charges
        ff_m=farfield(spec_s,sig); % Farfield for whole sphere        
        ff_ma=farfield(spec,sig); % Farfield for angular          
        %-------------------------% ellipticity for whole sphere
        for i=1:length(ff_m.e(1,1,:))
            fv=[mean(ff_m.e(:,1,i));mean(ff_m.e(:,2,i))];  
            [delta_circ_m(i,ien),n_circ_m(i,ien)]=polellip(fv);
        end        
        %-------------------------% ellipticity for angular
        for i=1:length(ff_ma.e(1,1,:))
            fv=[mean(ff_ma.e(:,1,i));mean(ff_ma.e(:,2,i))];
            [delta_circ_ma(i,ien),n_circ_ma(i,ien)]=polellip(fv);
            % Uncomment for phase
            %ExPhi = angle(mean(ff_ma.e(:,1,i)));
            %EyPhi = angle(mean(ff_ma.e(:,2,i)));
            %phi(i,ien) = EyPhi-ExPhi;   
        end              
        multiWaitbar( 'Ellipticity maps', ien / numel( enei2 ) );
    end  
    % Uncomment for phase
    %figure
    %prob_p = reshape( phi( :, ien ), size( xe ) ); %  Reshape CL spectra
    %imagesc( xx, yy, rad2deg(prob_p)  ); %  Add second part of structure
    %caxis([-90 90])
    multiWaitbar( 'CloseAll' ); %  Close waitbar
%%  Plot-section 2: CL maps 
    delete('data/ellipticity_map\*') % Delete old files  
    figure
    for ien = 1 : length(enei2)
        prob_n = reshape( n_circ_m( :, ien ), size( xe ) ); %  Reshape CL spectra
        prob_d = reshape( delta_circ_m( :, ien ), size( xe ) ); %  Reshape CL spectra
        subplot(2,2,1); %  Defines subplot window
        imagesc( xx, yy, prob_n  ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        colorbar
        if Double_en==0
                pbaspect([1 1 1])
            else
                pbaspect([2 1 1])
        end
        title(name(1));
        subplot(2,2,2); %  Defines subplot window
        imagesc( xx, yy, prob_d ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        colorbar
        if Double_en==0
                pbaspect([1 1 1])
            else
                pbaspect([2 1 1])
        end
        title(name(2));
        prob_na = reshape( n_circ_ma( :, ien ), size( xe ) ); %  Reshape CL spectra
        prob_da = reshape( delta_circ_ma( :, ien ), size( xe ) ); %  Reshape CL spectra
        subplot(2,2,3); %  Defines subplot window
        imagesc( xx, yy, prob_na ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        colorbar
        if Double_en==0
                pbaspect([1 1 1])
            else
                pbaspect([2 1 1])
        end
        title(name(3));
        subplot(2,2,4); %  Defines subplot window
        imagesc( xx, yy, prob_da ); %  Add second part of structure
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        colorbar
        if Double_en==0
                pbaspect([1 1 1])
            else
                pbaspect([2 1 1])
        end
        title(name(4));
        saveas(gcf,['data/ellipticity_map/' sprintf( '%04d', ien ) 'Ellipt_' num2str(ene2(ien)) 'eV' '.png'])
        saveas(gcf,['data/ellipticity_map/' sprintf( '%04d', ien ) 'Ellipt_' num2str(ene2(ien)) 'eV' '.fig'])
    end
    close
    imagefiles = dir('data/ellipticity_map/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/ellipticity_map/' currentfilename]));
        images{ii} = currentimage; % Load image files
    end
    f=figure;     
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
        name='Ellipticity ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end