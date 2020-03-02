%% Callback subfunctions to support UI actions for png single image
function slider1_callback(~,~,vars_b,vars_im,ene2, name)
    % Run slider1 which controls value of epsilon
    plotterfcn(vars_b,vars_im,ene2, name)
end