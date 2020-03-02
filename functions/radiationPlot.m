%% Radiation plot function, plots radiation results
function radiationPlot(dsca, ~, ~, dsca_s, dsca_sp, dsca_ss, ene, ien, detRad)
    %material dull 
    %shading interp
    %lighting none
    % ----- whole
    subplot( 2, 2, 1 );
    plot( trispherescale( dsca_s.p, dsca_s.dsca( :), 1 ), dsca_s.dsca( : ) );
    pbaspect([1 1 1])
    axis([-detRad detRad -detRad detRad -detRad detRad])
    xlabel( 'x (nm)' ); 
    ylabel( 'y (nm)' );
    zlabel( 'z (nm)' );
    grid on;
    material dull 
    shading interp
    lighting none
    colormap('jet')
    title(sprintf(['Radiation ', num2str(ene(ien)), ' eV ']));
    
    subplot( 2, 2, 2 );
    plot( trispherescale( dsca_sp.p, dsca_sp.dsca( : ), 1 ), dsca_sp.dsca( :) );
    pbaspect([1 1 1])
    axis([-detRad detRad -detRad detRad -detRad detRad])
    xlabel( 'x (nm)' ); 
    ylabel( 'y (nm)' );
    zlabel( 'z (nm)' );
    grid on;
    material dull 
    shading interp
    lighting none
    colormap('jet')
    title(sprintf(['Pol p ', num2str(ene(ien)), ' eV ']));
    
    subplot( 2, 2, 4 );
    plot( trispherescale( dsca_ss.p, dsca_ss.dsca( : ), 1 ), dsca_ss.dsca( : ) );
    pbaspect([1 1 1])
    axis([-detRad detRad -detRad detRad -detRad detRad])
    xlabel( 'x (nm)' ); 
    ylabel( 'y (nm)' );
    zlabel( 'z (nm)' );
    grid on;
    material dull 
    shading interp
    lighting none
    colormap('jet')
    title(sprintf(['Pol s ', num2str(ene(ien)), ' eV ']));
    
    % ----- Angular
    subplot( 2, 2, 3 );
    plot( dsca.p, dsca.dsca );
    pbaspect([1 1 1])
    axis([-detRad detRad -detRad detRad -detRad detRad])
    xlabel( 'x (nm)' ); 
    ylabel( 'y (nm)' );
    zlabel( 'z (nm)' );
    grid on;
    material dull 
    shading interp
    lighting none
    colormap('jet')
    title(sprintf(['Directional ', num2str(ene(ien)), ' eV ']));
    %subplot( 2, 3, 5 )
    %plot( trispherescale( dscap.p, dscap.dsca( : ), 1 ), dscap.dsca( : ) );
    %pbaspect([1 1 1])
    %axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
    %xlabel( 'x (nm)' ); 
    %ylabel( 'y (nm)' );
    %zlabel( 'z (nm)' );
    %grid on;
    %title(sprintf(['Angular p ', num2str(ene(ien)), ' eV ']));
    %subplot( 2, 3, 6 )
    %plot( trispherescale( dscas.p, dscas.dsca( : ), 1 ), dscas.dsca( : ) );
    %pbaspect([1 1 1])
    %axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
    %xlabel( 'x (nm)' ); 
    %ylabel( 'y (nm)' );
    %zlabel( 'z (nm)' );
    %grid on;
    %title(sprintf(['Angular s ', num2str(ene(ien)), ' eV ']));
    saveas(gcf,['data/radiation/' sprintf( '%04d', ien ) 'Radiation_' num2str(ene(ien)) 'eV' '.fig'])
    saveas(gcf,['data/radiation/' sprintf( '%04d', ien ) 'Radiation_' num2str(ene(ien)) 'eV' '.png'])
end
