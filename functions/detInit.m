%% Detector Init, Defines detector properties
function [spec,spec_s] = detInit(op,BEM_op,d_angle1,d_range1,d_angle2,detQ, detRad)
    %if d_angle2==0 %Special case for 0 degrees
        %angle_d1=[-pi, pi]; % Define detector sphere position 2
        %angle_d2=[0, d_range2/2]; % Define detector sphere position 2
    %elseif d_angle2==pi %Special case for 180 degrees
        %angle_d1=[-pi, pi]; % Define detector sphere position 2
        %angle_d2=[pi-d_range2/2, pi]; % Define detector sphere position 2
    %else
        %angle_d1=[d_angle1-d_range1/2, d_angle1+d_range1/2]; % Define detector sphere pos 1 Azimuth
        %angle_d2=[d_angle2-d_range2/2, d_angle2+d_range2/2]; % Define detector sphere position 2 polar
    %end
    angle_d1=[-pi, pi]; % Define detector sphere position 2
    angle_d2=[0, d_range1/2]; % Define detector sphere position 2
    if BEM_op==0 % (Retartded)
        pinfty_s = trisphere( detQ, detRad*2, 'interp', 'curv' ); % Sphere at selected range
        %pinfty_s = trispheresegment( linspace( 0, 2 * pi, 51 ), linspace( 0, pi, 51 ) );
        spec_s=spectrumret( pinfty_s, op ); % CL infinate sphere init
        pinfty_tmp = trispheresegment( linspace( angle_d1(1),  angle_d1(2), 51 ), linspace( angle_d2(1),  angle_d2(2),51 ), detRad*2, 'interp', 'curv'  ); % Sphere segment at infinity
        pinfty_tmprot=(rot(pinfty_tmp,radtodeg(d_angle2),[0,1,0]));
        pinfty=(rot(pinfty_tmprot,radtodeg(d_angle1),[0,0,1]));
        spec = spectrumret( pinfty, op ); % Set up spectrum object
    else % (Quasistatic)        
        pinfty_s = trisphere( detQ, detRad*2, 'interp', 'curv' ); % Sphere at selected range
        spec_s=spectrum( pinfty_s, op ); % CL infinate sphere init
        pinfty_tmp = trispheresegment( linspace( angle_d1(1),  angle_d1(2), 51 ), linspace( angle_d2(1),  angle_d2(2), 51 ), detRad*2, 'interp', 'curv'  ); % Sphere segment at infinity
        pinfty_tmprot=(rot(pinfty_tmp,radtodeg(d_angle2),[0,1,0]));
        pinfty=(rot(pinfty_tmprot,radtodeg(d_angle1),[0,0,1]));
        spec = spectrum( pinfty, op ); % Set up spectrum object
    end  
end