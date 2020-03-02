%% BEM solver, computes electron beam excitation
function [bem, exc, exc_CL, exc_CL_M, exc_EELS_M]= BEMSolver(p, op, pw_pol, imp,impact,imp_spec, width, vel, Source_op, dip_p, dip_d, d_range2a)
    bem = bemsolver( p, op ); %  BEM solver
    if Source_op==1 % Excitation source selection (1 for e-beam, 2 for laser, 3 for dipole)
        exc = electronbeam( p, imp( : ) * [ 1, 0 ], width, vel, op ); %  Electron beam excitation
        exc_CL = electronbeam( p, imp_spec, width, vel, op );%  Electron beam excitation for specific imp
        exc_CL_M = electronbeam( p, impact, width, vel, op ); %  CL map
        exc_EELS_M = electronbeam( p, impact, width, vel, op ); %  EELS excitation map
    elseif Source_op==2 % Laser    
        if pw_pol==1
            exc_CL=planewave( [ 1, 0, 0 ], [ 0, 0, 1], op ); %  p pol
        elseif pw_pol==2
            exc_CL=planewave( [ 0, 1, 0 ], [ 0, 0, 1], op ); %  s pol
        elseif pw_pol==3
            exc_CL=planewave( [ 0.5, 0.5, 0 ], [ 0, 0, 1], op ); %  s pol
        elseif pw_pol==4
            if d_range2a<-45
                x_pol=(90-abs(d_range2a))/45/2;
                y_pol=-1+x_pol;
            elseif d_range2a<0
                y_pol=d_range2a/45/2;
                x_pol=1-abs(y_pol);                
            elseif abs(d_range2a)<=45                
                y_pol=d_range2a/45/2;
                x_pol=1-y_pol;
            elseif abs(d_range2a)>45
                x_pol=(90-abs(d_range2a))/45/2;
                y_pol=1-x_pol;     
            end
            exc_CL=planewave( [ x_pol, y_pol, 0 ], [ 0, 0, 1], op ); %  rotation pol
        end
        %------------- Polarization rotation
        j=-90;
        for i=1:1:37
            if j<-45
                x_pol=(90-abs(j))/45/2;
                y_pol=-1+x_pol;
            elseif j<0
                y_pol=j/45/2;
                x_pol=1-abs(y_pol);                
            elseif abs(j)<=45                
                y_pol=j/45/2;
                x_pol=1-y_pol;
            elseif abs(j)>45
                x_pol=(90-abs(j))/45/2;
                y_pol=1-x_pol;     
            end
            j=j+5;
            pol_m(i,:)=[x_pol, y_pol, 0];
        end
        exc = planewave( pol_m, [ 0, 0, 1], op );
        exc_CL_M=0;
        exc_EELS_M=0;
    else % Dipole 
        %------------- Single position
        pt = compoint( p, [ dip_p(1), dip_p(2), dip_p(3)] ); %  Dipole position
        if dip_d==1
            exc_CL = dipole( pt, [ 1, 0, 0], op ); %  Dipole excitation X
        elseif dip_d==2
            exc_CL = dipole( pt, [ 0, 1, 0], op ); %  Dipole excitation Y
        else
            exc_CL = dipole( pt, [ 0, 0, 1], op ); %  Dipole excitation Z
        end    
        %------------- for multiple z values
        mult=ones(length(imp),1); 
        pt = compoint( p, [ imp( : ), dip_p(2)*mult, dip_p(3)*mult] ); %  Dipole position
        if dip_d==1
            exc = dipole( pt, [ 1, 0, 0], op ); %  Dipole excitation X
        elseif dip_d==2
            exc = dipole( pt, [ 0, 1, 0], op ); %  Dipole excitation Y
        else
            exc = dipole( pt, [ 0, 0, 1], op ); %  Dipole excitation Z
        end
        %------------- for mapping
        mult=ones(length(impact(:,1)),1); 
        pt = compoint( p, [ impact(:,1), impact(:,2), dip_p(3).*mult ]); %  Dipole position
        if dip_d==1
            exc_CL_M = dipole( pt, [ 1, 0, 0], op ); %  Dipole excitation X
        elseif dip_d==2
            exc_CL_M = dipole( pt, [ 0, 1, 0], op ); %  Dipole excitation Y
        else
            exc_CL_M = dipole( pt, [ 0, 0, 1], op ); %  Dipole excitation Z
        end
        %------------- for multiple z values (3D direction)
        mult=ones(length(imp),1); % for multiple z values
        pt = compoint( p, [ imp( : ), dip_p(2)*mult, dip_p(3)*mult]); %  Dipole position
        exc_EELS_M=dipole( pt, op ); %  Dipole excitation 3D
    end
end