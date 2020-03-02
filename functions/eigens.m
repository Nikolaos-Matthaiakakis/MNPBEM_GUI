%%  Plasmonic Eigenmodes
function eigens(BEM_op, p, neg, op, bem)
    if BEM_op==1
        %[ ene_pe, ur ] = plasmonmode( p, neg, op );
        %ene_pe_p=abs(ene_pe);
        %enei_pe = eV2nm ./ ene_pe_p;
        figure % Plot eigenmodes
        plot( p, bem.ur);
        xlabel( 'x (nm)' ); 
        ylabel( 'y (nm)' );
        zlabel( 'z (nm)' );
        title( 'Eigenmodes' );
        material dull 
        shading interp
        lighting none
        %colormap (jet)
        saveas(gcf,'data/Eigen.fig')
        saveas(gcf,'data/Eigen.png')
    end
end