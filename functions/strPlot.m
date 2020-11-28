%% Nanostructure for GUI
function strPlot()
    % Load input variables  
    mat = dir('tmp/*.mat'); 
    for q = 1:length(mat) 
        load(sprintf(['tmp/' mat(q).name])); 
    end 
    dip_p=[dip_x dip_y dip_z ];
    neg=1;

    %% Structure input
    ax_el = [ d_length, d_length2, d_length3;  d_length, d_length2, d_length3;  d_length, d_length2, d_length3 ]; %  Axes of ellipsoid nm

    %% Refractive index Init
    epstab = refInit(Double_en,Cover_en, eps_env, eps_part, eps_str, eps_part2, eps_p, eps_list_en, eps_list_en2, substrate, substrate2, eps_p2, eps_part3, eps_list_en3, Source_op, qd_en);

    %% BEM initialization
    op = BEMInit(BEM_op,neg, intPoint, dIntPoint, Relcut);
   
    %% Selection of nanostructure
    if typeStr==11 % User defined
        p_loc= strInitUser(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en);    
    else
        p_loc = strInit(typeStr, Double_en, d_length, edge_profile, Cover_en, epstab, op, ax_el, d_c, shift_l, d_length2, d_length3, angles_n, grid_1, grid_2, hdata, substrate, sub_t, substrate2, sub_t2, xtilt, ytilt, ztilt, dip_p, Source_op, qd_en);
    end
    f_tmp = figure('visible','off');
    %p_loc = flip( p_loc, 3 ); 
    plot(p_loc, 'FaceColor', 0.8 * [ 1, 1, 1 ], 'EdgeColor','b')
    hold on
    pos1 = quad(p_loc);
    plot3( pos1( :, 1 ), pos1( :, 2 ), pos1( :, 3 ), 'r.' ,'MarkerSize',3 );
    axis tight
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    saveas(f_tmp,'data/str_tmp.png')
    % 3D plot
    f_tmp2 = figure('visible','off');
    plot(p_loc, 'FaceColor', 0.8 * [ 1, 1, 1 ], 'EdgeColor','b')
    hold on
    pos1 = quad(p_loc);
    plot3( pos1( :, 1 ), pos1( :, 2 ), pos1( :, 3 ), 'r.' ,'MarkerSize',1 );
    axis on
    pbaspect([1 1 1])
    %if Double_en==0
            %pbaspect([1 1 1])
            %axis([-d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5 -d_length/1.5-d_length/5 d_length/1.5+d_length/5])
        %else
            %pbaspect([2 1 1])
            %axis_tmp=(2*d_length+shift_l)/1.6;
            %axis([-axis_tmp axis_tmp -axis_tmp/2 axis_tmp/2 -axis_tmp/2 axis_tmp/2])
    %end
    title('Structure geometry')
    material dull 
    shading interp
    lighting none
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    saveas(f_tmp2,'data/str.fig')
end