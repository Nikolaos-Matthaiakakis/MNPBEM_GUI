%% Beam Init, loads e-beam properties for the simulation
function [width, vel, enei2, impact, xx, yy, xe, ye] = beamInit( width_d, e_beam, ene2, Double_en, d_length, res_map, XYField, shift_l)
    [ width, vel ] = deal( width_d, eelsbase.ene2vel( e_beam ) ); %  Width of electron beam and electron velocity
    units;  
    enei2 = eV2nm ./ ene2;

    %% Maps settings
    if XYField ~= 0
        d_length=XYField;
    end
    if Double_en==0
        %[ xe, ye ] = meshgrid( linspace( -d_length/1.5-d_length/5, d_length/1.5+d_length/5, res_map ), linspace( 0, d_length/1.5+d_length/5, res_map/2 ) ); %  Mesh for electron beams
    [ xe, ye ] = meshgrid( linspace( -d_length/1.5-d_length/5, d_length/1.5+d_length/5, res_map ), linspace( -d_length/1.5-d_length/5, d_length/1.5+d_length/5, res_map/2 ) ); %  Mesh for electron beams
    else
        beam_temp=(2*d_length+shift_l)/1.6;
        [ xe, ye ] = meshgrid( linspace( -beam_temp, beam_temp, res_map ), linspace( -beam_temp/2 , beam_temp/2, res_map/2 ) ); %  Mesh for electron beams
    end
    impact = [ xe( : ), ye( : ) ]; %  Impact parameters
    xx = [   min( xe( : ) ), max( xe( : ) ) ]; %  x and y limits
    %yy = [ - max( ye( : ) ), max( ye( : ) ) ]; %  x and y limits 
    yy = [   min( ye( : ) ), max( ye( : ) ) ]; %  x and y limits
end  

