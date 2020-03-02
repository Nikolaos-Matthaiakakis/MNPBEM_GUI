function Decay_rad(op, epstab, d_length, typeStr, exc_CL, bem, ene2, enei2, dip_p, dip_d, p)
    multiWaitbar( 'Total and radiative decay', 0, 'Color', 'g', 'CanCancel', 'on' );
    if typeStr==0 %  Solution for sphere
        mie = miesolver( epstab{ 2 }, epstab{ 1 }, d_length, op, 'lmax', 40 ); %  Mie solver
    end
    for ien=1:length(enei2)
        cancelRun() % Cancel function
        sig = bem \ exc_CL(p, enei2( ien ) ); %  Surface charges 
        [ tot( ien, : ), rad( ien, : ) ] = exc_CL.decayrate( sig ); %  Total and radiative decay rate
        if typeStr==0 %  Solution for sphere
            [ tot0( ien, : ), rad0( ien, : ) ] = mie.decayrate( enei2( ien ),dip_p(3) );
            if dip_d==3
                tot0p( ien, : )=tot0( ien, 2 );
                rad0p( ien, : )=rad0( ien, 2 );
            else
                tot0p( ien, : )=tot0( ien, 1 );
                rad0p( ien, : )=rad0( ien, 1 );
            end
        end
        multiWaitbar( 'Total and radiative decay', ien / numel( enei2 ) );
    end
    
    multiWaitbar( 'CloseAll' );
    %%  final plot
    figure
    subplot(2,1,1)
    if typeStr==0 %  Solution for sphere
        plot( ene2, tot, ene2, tot0p, 'g*'); 
        legend('Total-simulation','Total-Mie')
    else
        plot( ene2, tot)
    end    
    xlabel( 'Energy (eV)' );
    ylabel( 'Total decay rate' );
    subplot(2,1,2)
    if typeStr==0 %  Solution for sphere
        plot( ene2, rad, ene2, rad0p, 'g*'); 
        legend('Total-simulation','Total-Mie')
    else
        plot( ene2, rad)
    end    
    xlabel( 'Energy (eV)' );
    ylabel( 'Radiative decay rate' );
    saveas(gcf,'data/decay.fig')
    saveas(gcf,'data/decay.png')
end
