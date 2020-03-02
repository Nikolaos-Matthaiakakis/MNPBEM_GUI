%% Callback plotter for single figure
function plotterfcnFig(vars_b, f, dire)
    delete(findall(findall(f,'Type','axe'),'Type','text')) % Deletes previous axis values
    cla 
    delete(gca) % Clears previous graphical element 
    imagefiles = dir(sprintf(['data/' dire '/*.fig']));  % Find files        
    g = openfig(sprintf(['data/' dire '/' imagefiles(round(get(vars_b.b,'Value'))).name]), 'invisible'); % Loads new graphic
    copyobj(get(g, 'CurrentAxes'), f); % Refreshes figure with new axis and graphical element
    delete(g);  
end