%% Callback plotter for single png figure
function plotterfcn(vars_b,vars_im,ene2, name)    
    imshow(vars_im(round(get(vars_b.b,'Value'))).images) % Plots new image
    title([name num2str(ene2(round(get(vars_b.b,'Value')))) 'eV'] );
end