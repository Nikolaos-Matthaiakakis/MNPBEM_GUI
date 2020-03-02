%% Callback plotter for subplot figure
function plotterfcnFig3(vars_b, f, dire)
    % Plots the image
    fx=f.Children;
    if length(fx)>=5
        for i=1:4 % Deletes all subplot axis and graphical elements
            delete(findall(findall(f,'Type','axe'),'Type','text'))
            if length(fx)>=5
                delete(f.Children(length(fx)-3))
            end
        end
        imagefiles = dir(sprintf(['data/' dire '/*.fig']));  % Find files 
        g = openfig(sprintf(['data/' dire '/' imagefiles(round(get(vars_b.b,'Value'))).name]), 'invisible'); 
        subFig2 = findobj( g, 'Type', 'Axes' ); % Finds subplot indexes
        copyobj(subFig2, f); % Plots new image axis and graphical elements
        delete(g);
    end
end