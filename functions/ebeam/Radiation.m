%%  Radiation calculation   
function Radiation( exc_CL, bem, spec, spec_s,  ene, enei, detRad)
    multiWaitbar( 'Radiation', 0, 'Color', 'g', 'CanCancel', 'on' ); %  Loop over wavelengths
    delete('data/radiation\*') % Delete old files
    dire='radiation';
    for ien = 1 : length( enei )
        drawnow() % Reads new inputs
        cancelRun() % Cancel function
        sig = bem \ exc_CL( enei( ien ) ); %  Surface charges
        
        % ----- Whole
        ff_s=farfield(spec_s,sig); % Farfield for whole sphere
        [~,dsca_s] = scattering(ff_s); %  Farfield scattering for whole sphere
        %-------------------------% p pol
        ffp_s=ff_s;
        ffp_s.e(:,2)=0; % E field in x direction set to 0
        ffp_s.h(:,1)=0; % B field in y direction set to 0
        [~,dsca_sp]=scattering(ffp_s); %  Farfield calculation
        %-------------------------% s pol
        ffs_s=ff_s;
        ffs_s.e(:,1)=0; % E field in y direction set to 0
        ffs_s.h(:,2)=0; % B field in x direction set to 0
        [~,dsca_ss]=scattering(ffs_s); %  Farfield calculation        
        
        % ----- Angular
        ff=farfield(spec,sig); % Farfield angular
        [~,dsca] = scattering(ff); %  Farfield scattering
        %-------------------------% p pol
        ffp=ff;
        ffp.e(:,2)=0; % E field in x direction set to 0
        ffp.h(:,1)=0; % B field in y direction set to 0
        [~,dscap]=scattering(ffp); %  Farfield calculation
        %-------------------------% s pol
        ffs=ff;
        ffs.e(:,1)=0; % E field in y direction set to 0
        ffs.h(:,2)=0; % B field in x direction set to 0
        [~,dscas]=scattering(ffs); %  Farfield calculation
        f=figure;
        radiationPlot(dsca, dscap, dscas, dsca_s, dsca_sp, dsca_ss, ene, ien, detRad)
        close(f)
        multiWaitbar( 'Radiation', ien / numel( enei ) );
    end    
    multiWaitbar( 'CloseAll' ); %  Close waitbar
    imagefiles = dir('data/radiation/*.fig');  % Find files 
    f=openfig(sprintf(['data/radiation/' imagefiles(1).name]), 'visible');
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
        set(be,'Callback',{@slider1_callbackFig3, vars_b, f, dire});
        plotterfcnFig3(vars_b, f, dire)         
    end
end