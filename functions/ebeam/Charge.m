%%  Charge calculation   
function Charge( exc_CL, bem, Double_en, d_length,  ene, enei, BEM_op)
    multiWaitbar( 'Charge', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths
    delete('data/charge\*') % Delete old files
    dire='charge';    
    f=figure;
    for ien = 1 : length( enei )
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL( enei( ien ) ); %  Surface charges
        grid on
        if BEM_op==0 || BEM_op==2
            plot (sig.p,sig.sig1); %  Plots structure, mesh, and charges
        else
            plot (sig.p,sig.sig); %  Plots structure, mesh, and charges
        end
        material dull 
        shading interp
        lighting none
        colormap (jet)        
        if Double_en==0
            pbaspect([1 1 1])
            axis([-d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5])
        else
            pbaspect([2 1 1])
            axis_tmp=d_length+d_length/1.5+d_length/2.5+d_length/1.5+d_length/5;
            axis([-d_length-d_length/1.5-d_length/2.5 d_length/1.5+d_length/5 -axis_tmp/4 axis_tmp/4 -axis_tmp/4 axis_tmp/4])
        end
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        zlabel( 'z (nm)' );
        title(sprintf([ 'Charge for ' num2str(ene(ien)) ' eV' ]));
        saveas(gcf, ['data/charge/' sprintf( '%04d', ien ) 'Charge_' num2str(ene(ien)) 'eV' '.fig'])
        saveas(gcf, ['data/charge/' sprintf( '%04d', ien ) 'Charge_' num2str(ene(ien)) 'eV' '.png'])
        multiWaitbar( 'Charge', ien / numel( enei ) );
    end
    multiWaitbar( 'CloseAll' ); %  Close waitbar
    val=1;
    % Assign slider properties
    if length(ene)>1
        be = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'Value',val, 'Min',1, 'Max',length(ene), 'SliderStep', [1/(length(ene)-1) 1/(length(ene)-1)]);
        bgecolor = f.Color;
        bel1 = uicontrol('Parent',f,'Style','text','Position',[50,24,23,23],...
                'String',min(ene),'BackgroundColor',bgecolor);
        bel2 = uicontrol('Parent',f,'Style','text','Position',[500,24,23,23],...
                'String',max(ene),'BackgroundColor',bgecolor);
        bel3 = uicontrol('Parent',f,'Style','text','Position',[240,5,100,13],...
                'String','eV','BackgroundColor',bgecolor);    
        % Set up callbacks
        vars_b=struct('b',be);
        set(be,'Callback',{@slider1_callbackFig, vars_b, f, dire});
        plotterfcnFig(vars_b, f, dire)     
    end
end