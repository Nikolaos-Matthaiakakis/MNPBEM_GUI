%%  Field calculation, reruns electronbeam, uses enei2 instead of enei
function field_cal(d_length, ene2, enei2, Double_en, bem, exc_CL, BEM_op, p, op, res_map, XYField, cmSt, cmEn, shift_l) 
    delete('data/field\*') % Delete old files
    if XYField ~= 0
        d_length=XYField;
    end
    if Double_en==0    
        [ x, z ] = meshgrid( linspace( - d_length/1.5-d_length/5, d_length/1.5+d_length/5, res_map ) );  %  Mesh for calculation of electric field XY
    else
        axis_tmp=(2*d_length+shift_l)/1.6;
        [ x, z ] = meshgrid( linspace( -axis_tmp, axis_tmp, res_map ), linspace(  -axis_tmp/2, axis_tmp/2, res_map/2 ) ); %  Mesh for electron beams
    end
    %    object for electric field
    %    MINDIST controls the minimal distance of the field points to the
    %    particle boundary     
    multiWaitbar( 'Field calculation', 0, 'Color', 'g', 'CanCancel', 'on' );
    for ien = 1 : length( enei2 ) %  Loop over energies and plots charge
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL( enei2( ien ) ); 
        %-------------------------%  Calculate for XY plane
        cancelRun()
        emesh = meshfield( p, x, z, 0 , op, 'mindist', 0.1 );
        e = emesh( sig ); %  Induced and incoming electric field  
        %e = emesh( sig )+ emesh( exc_CL.field( emesh.pt, enei2(ien) ) );  %  Induced and incoming electric field
        ee = sqrt( dot( e, e, 3 ) ); %  Norm of electric field
        %-------------------------% Calculate for XZ plane
        cancelRun()
        emesh2 = meshfield( p, x, 0, z , op, 'mindist', 0.1 ); %  Mesh for calculation of electric field XZ    
        e2 = emesh2( sig ); %  Induced and incoming electric field  
        %e2 = emesh2( sig )+ emesh2( exc_CL.field( emesh2.pt, enei2(ien) ) ); %  Induced and incoming electric field
        ee2 = sqrt( dot( e2, e2, 3 ) ); %  Norm of electric field          
        %-------------------------% Plot figure
        figure('Position', [100, 100, 1000, 500]);        
        %-------------------------% Plot for XY plane
        ax2=subplot( 1, 2, 1 );
        imagesc( x( : ), z( : ), log10(ee));
        hold on;  
        plot(ax2, exc_CL.impact(1) , exc_CL.impact(2) , '*w:' ); % plot electron trajectory
        hold off;
        if Double_en==0 %  For single structure
            pbaspect([1 1 1])
            axis([-inf inf -inf inf])
        else %  For double structure
            pbaspect([2 1 1])
            axis([-inf inf -axis_tmp/2 axis_tmp/2])
        end
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        if cmSt ~= 0 || cmEn ~= 0
            caxis([cmSt cmEn]); 
        else
            caxis('auto')
        end
        colorbar; colormap (ax2,hot( 255 ));
        title(sprintf([ 'Field XY ' num2str(ene2(ien)) ' eV (log)' ]));
        %-------------------------% Plot for XZ plane  
        ax3=subplot( 1, 2, 2 );
        imagesc( x( : ), z( : ), log10(ee2) );
        hold on;  
        plot(ax3, exc_CL.impact(1) * [ 1, 1 ], (d_length/1.5+d_length/5) * [ - 1, 1 ], 'w:' ); % plot electron trajectory
        hold off;
        if Double_en==0
            pbaspect([1 1 1])
            axis([-inf inf -inf inf])
        else
            pbaspect([2 1 1])
            axis([-inf inf -axis_tmp/2 axis_tmp/2])
        end
        set( gca, 'YDir', 'norm' );
        xlabel( 'x (nm)' ); 
        ylabel( 'z (nm)' );
        if cmSt ~= 0 || cmEn ~= 0
            caxis([cmSt cmEn]); 
        else
            caxis('auto')
        end
        colorbar; colormap (ax3,hot( 255 ));
        title(sprintf([ 'Field XZ ' num2str(ene2(ien)) ' eV (log)' ]));
        multiWaitbar( 'Field calculation', ien / numel( enei2 ) );
        saveas(gcf,sprintf(['data/field/' sprintf( '%04d', ien ) 'Field_' num2str(ene2(ien)) 'eV' '.fig']))
        saveas(gcf,sprintf(['data/field/' sprintf( '%04d', ien ) 'Field_' num2str(ene2(ien)) 'eV' '.png']))      
        close
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar    
    imagefiles = dir('data/field/*.png') ; % Find files in folder       
    nfiles = length(imagefiles); % Number of files found
    for ii=1:nfiles
        currentfilename = imagefiles(ii).name;
        currentimage = imread(sprintf(['data/field/' currentfilename]));
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
        name='Field ';
        set(b,'Callback',{@slider1_callback,vars_b,vars_im,ene2, name});
        plotterfcn(vars_b,vars_im,ene2, name)  
    else
        imshow(images{1})
    end    
end