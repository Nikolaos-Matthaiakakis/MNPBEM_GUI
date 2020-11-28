%  -In order for the code to work it has to be placed in the function folder
%  environment.  
%  http://physik.uni-graz.at/mnpbem/download.php
%-------------------------% 
%  -If this code is used in a publication please include the references as
%  mentioned in the following website http://physik.uni-graz.at/mnpbem/#2
%-------------------------% 
%  -For more information please visit http://physik.uni-graz.at/mnpbem/
%-------------------------% 

%-------------------------% 
%  -Main program section does not need to be changed other than for specific
%  reasons.
%-------------------------% 
% For user defined structures please read the strInitUser_readme file in the function/user folder
%-------------------------% 
% For user defined permittivity please read the userNK_readme file in the @epstable folder
%-------------------------%
% For user defined function please read the user_fun file in the
% function/user folder and edit line 265 user_fun() in this m file
%-------------------------%

%  -Notes:
%  plot(p,'EdgeColor','b','nvec',1). % Run to verify face orientation
%  prad(:,ien) = exc.rad( sig ); 
%  plot(ene,  mie.loss( imp-d_length/2, enei, vel ),'*-'); %plots Mie loss
%  plot( ene, mie.sca( enei ), '*-' ); %plots Mie
%  [ eps, k ] = epsobj( enei ); evaluation of dielectric function at wavelength ENEI


%  In the BEM approach, outside the sphere the surface charges and currents sig2 and h2 account
%  for the induced fields (we have to add the external fields to obtain the total fields), 
%  while inside the particle sig1 and h1 account for the total fields.

%  for polarization compute farfields and then set, for instance, the x-component of the electric field
%  (and the corresponding components of the magnetic field) equal to
%  zero. Then continue as normal.

%  Electron beam size should be smaller than average face length

%  Electron beam propagates in the positive direction

% Author Nikolaos Matthaiakakis CL_GUI v2
function CL_app()
tic % Timer on
    % Change the current folder to the folder of this m-file.
    % Load input variables   
    mat = dir('tmp/*.mat');
    for q = 1:length(mat) 
        load(sprintf(['tmp/' mat(q).name])); 
    end   

    %% Default colormap
    set(0,'DefaultFigureColormap',feval('jet')); % Selects colormap for 3D plots

    %% General solver and analysis settings
    neg=30; % Number of Eigenmodes for plasmonic Eigenmode calculation

    %% Structure input
    ax_el = [ d_length, d_length2, d_length3;  d_length, d_length2, d_length3;  d_length, d_length2, d_length3 ]; %  Axes of ellipsoid nm
     
    %% Detector input
    d_angle1=degtorad(d_angle1a); % XY plane angle
    d_angle2=degtorad(180-d_angle2a); % XZ angle
    d_range1=degtorad(d_range1a); % Vertical to Z angle range
    
    %% Excitation source input
    ene2= linspace( minene2, maxene2, stepene2 );%2 2.162 2.53 2.8 3.95] ; %  Energies for CL map and EELS map sections 
    %-------------------------% 
    imp=linspace( minimp, maximp, stepimp );%,linspace( minimpY, maximpY, stepimp ) ]; %  Beam position or range of positions for EELS loss and CL spectra nm (Mie does not work for imp < diameter value)
    imp_spec=[imp_specX, imp_specY]; % XY beam values for specific position
    %-------------------------% 
    dip_p=[dip_x dip_y dip_z ];
    %-------------------------% 
    d_range2=degtorad(d_range2a); % light polarization
    
    %% Main program

    %% Init
    
    %% Refractive index Init
    epstab = refInit(Double_en,Cover_en, eps_env, eps_part, eps_str, eps_part2, eps_p, eps_list_en, eps_list_en2, substrate, substrate2, eps_p2, eps_part3, eps_list_en3, Source_op, qd_en);

    %% BEM initialization
    op = BEMInit(BEM_op, neg, intPoint, dIntPoint, Relcut);   
    
    %% Beam Init
    [width, vel, enei2, impact, xx, yy, xe, ~] = beamInit( width_d, e_beam, ene2, Double_en, d_length, res_map, XYField, shift_l);

    %% Selection of nanostructure
    if typeStr==11 % User defined  
        if substrate==1
            [p, op, exc_flag, bem_tmp] = strInitUser_Al_comp(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en, ene2, Field_en, xe, xe);
        else
            [p, op] = strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en);
        end
        epstab=p.eps;
    else % GUI
        p = strInit(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en);
    end

    %% Detector Init
    [spec,spec_s] = detInit(op,BEM_op,d_angle1,d_range1,d_angle2,detQ,detRad);
    
    %%  BEM solution 
    if exist('exc_flag','var')% Check for user defined excitation      
        exc=exc_flag;
        exc_CL=exc_flag;
        bem=bem_tmp;
    else
        [bem, exc, exc_CL, exc_CL_M, exc_EELS_M]= BEMSolver(p, op, pw_pol, imp,impact,imp_spec, width, vel, Source_op, dip_p, dip_d, d_range2a);
    end
        %%  Plasmonic Eigenmodes
    if BEM_op==1
        eigens(BEM_op, p, neg, op, bem)
    end       
    
    %%  E-beam functions
    if Source_op==1
        %%  EELS loss and CL spectra calculation  
        if EELS_CL==1
            EELS_CL_calc(op, vel, exc, bem, Double_en, d_length, epstab, imp,  imp_spec, ene2, enei2, p, spec, spec_s, exc_CL, typeStr, shift_l)
        end

        %%  CL Angle Scan 
        if CLang==1
            CL_Angle(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, d_angle1, vel)
        end

        %%  Charge calculation  
        if charge_f==1
            Charge( exc_CL, bem, Double_en, d_length,  ene2, enei2, BEM_op)
        end

        %%  Radiation calculation  
        if rad_f==1
            Radiation( exc_CL, bem, spec, spec_s,  ene2, enei2, detRad)
        end

        %%  CL maps calculation
        if CL_map_en==1
            CL_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, vel)
        end

        %%  Ellipticity maps calculation
        if Pol_map_en==1
            Pol_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en)
        end

        %%  Ellipticity Angle Scan 
        if Pol_ang_en==1
            Pol_Angle(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, d_angle1)
        end

        %%  Field calculation, reruns electronbeam
        if Field_en==1 
            field_cal(d_length, ene2, enei2, Double_en, bem, exc_CL, BEM_op, p, op, res_map, XYField, cmSt, cmEn, shift_l);
        end

        %%  EELS maps calculation
        if EELS_en==1
            EELSmap( impact, ene2, enei2, exc_EELS_M, xx, yy, bem, xe, Double_en)
        end
        
        %%  Tilt scan calculation
        if tilt_en==1
            tilt_scan(op, ene2, enei2, p, imp_spec, width, vel, spec_s)
        end
        
    %%  Laser functions 
    elseif Source_op==2 
        %%  Scattering and farfield spectra calculation  
        if EELS_CL==1
            FF_sca_calc(op, exc_CL, bem, Double_en, d_length, epstab, ene2, enei2, p, spec, spec_s, typeStr, shift_l)
        end
        
        %%  Farfield Angle Scan 
        if CLang==1
            FF_Angle(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, p, d_angle1)
        end

        %%  Charge calculation  
        if charge_f==1
            Charge_l( exc_CL, bem, Double_en, d_length,  ene2, enei2, BEM_op, p)
        end

        %%  Radiation calculation  
        if rad_f==1
            Radiation_l( exc_CL, bem, spec, spec_s,  ene2, enei2, detRad, p)
        end

        %%  Ellipticity Angle Scan 
        if Pol_ang_en==1
            Pol_Angle_l(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, p, d_angle1)
        end

        %%  Field calculation
        if Field_en==1 
            field_cal_l(d_length, ene2, enei2, Double_en, bem, exc_CL, BEM_op, p, op, res_map, XYField, cmSt, cmEn, shift_l);
        end
        
        %%  Tilt scan calculation
        if tilt_en==1
            %tilt_scan_l(op, ene2, enei2, p, spec_s, pw_pol)
            tilt_scan_l(op, ene2, enei2, p, spec, pw_pol)
        end 
        
        %%  Polarization rotation calculation
        if pol_rot==1
        	pol_rot_scan(p, ene2, enei2, spec_s, exc, bem)
        end
               
    %%  Dipole functions     
    else
        %%  Scattering and farfield spectra calculation  
        if EELS_CL==1
            FF_sca_calc_d(op, exc_CL, bem, Double_en, d_length, epstab, ene2, enei2, p, spec, spec_s, typeStr, shift_l)
        end
        
        %%  Farfield Angle Scan 
        if CLang==1
            FF_Angle(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, p, d_angle1)
        end

        %%  Charge calculation  
        if charge_f==1
            Charge_l( exc_CL, bem, Double_en, d_length,  ene2, enei2, BEM_op, p)
        end

        %%  Radiation calculation  
        if rad_f==1
            Radiation_l( exc_CL, bem, spec, spec_s,  ene2, enei2, detRad, p)
        end

        %%  Ellipticity Angle Scan 
        if Pol_ang_en==1
            Pol_Angle_l(op, exc_CL, bem, ene2, enei2, BEM_op, d_range1, detQ, detRad, p, d_angle1)
        end

        %%  Field calculation
        if Field_en==1 
            field_cal_l(d_length, ene2, enei2, Double_en, bem, exc_CL, BEM_op, p, op, res_map, XYField, cmSt, cmEn, shift_l);
        end        
        
        %%  Decay rate calculation
        if decay_en==1
            Decay_rad(op, epstab, d_length, typeStr, exc_CL, bem, ene2, enei2, dip_p, dip_d, p)
        end  
        
        %%  Decay scan/LDOS 1D calculation 
        if LDOS_en==1 && qd_en==0
            Decay_scan( op, epstab, d_length, typeStr, exc, bem, ene2, enei2, dip_p, dip_d, p, imp)
        end
        
        %%  LDOS 3D calculation 
        if LDOS3D_en==1 && qd_en==0
            Decay_scan3D(  d_length, exc_EELS_M, bem, ene2, enei2, p, imp)
        end
        
        %%  Tilt scan calculation
        if tilt_en==1
            tilt_scan_d(op, ene2, enei2, p, spec_s, dip_d, dip_p)
        end 
        
        %%  Decay map calculation 
        if Dmaps_en==1 && qd_en==0
            Decay_map( exc_CL_M, bem, ene2, enei2, p, Double_en, xe, xx, yy)
        end
        
        %%  Luminescence maps calculation 
        if Ld_map_en==1 && qd_en==0
            L_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, p)
        end
        
        %%  Elipticity maps calculation 
        if Pold_map_en==1 && qd_en==0
            Pol_map_d(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, p)
        end
        
        %%  Luminescence enchancement spectra
        if Enh_sp_en==1
            if qd_en==0
                FF_sca_calc_den(exc_CL, bem, Double_en, d_length, ene2, enei2, p, spec, spec_s, shift_l)
            else
                FF_sca_qd_den(exc_CL, bem, Double_en, d_length, ene2, enei2, p, spec, spec_s, shift_l, op, dip_d)
            end
        end
        
        %%  Luminescence enhancement maps calculation
        if Enh_map_en==1 && qd_en==0
            Len_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en, p )
        end
    end
    
    %%  Area for user defined function
    %   Please enter your call function here as user_fun(input1,input2,...)
    %   and also in CL_GUI/functions/user
    if user_en==1
    	%user_fun(exc_CL, bem, Double_en, d_length,  ene2, enei2, BEM_op)
        
        %Chiral_CL_calc(op, vel, exc, bem, Double_en, d_length, epstab, imp,  imp_spec, ene2, enei2, p, spec, spec_s, exc_CL, typeStr)
        Chiral_map(impact, ene2, enei2, spec_s, spec, xe, bem, exc_CL_M, xx, yy, Double_en)
    end
    %%  End 
    toc %  End timer
end