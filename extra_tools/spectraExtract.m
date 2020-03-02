% Put in folder with CL maps to extract from
% Verify only map .fig files are in the folder, after running REMOVE
% produced .fig image to rerun !!!!!

clear all
clc
xc=0; %input X value nm
yc=160; %input Y value nm

name=["CL", "CL p", "CL s", "CL_a", "CL_a p", "CL_a s"];
imagefiles = dir(sprintf('*.fig'));  % Find files   

for i=1:length(imagefiles) % Runs for each .fig file
    evRangeStr=imagefiles(i).name(8:end-6); 
    evRange(i)=str2num(evRangeStr); % Loads eV range
    fig=openfig(imagefiles(i).name, 'invisible'); % Loads image
    subFig = findobj( gcf, 'Type', 'Axes' ); % Finds subplot indexes
    for j=1:6 % Runs for each subplot
        xAxes=-subFig(j).XLim(1)+subFig(j).XLim(2); % Finds total x range from axis
        yAxes=-subFig(j).YLim(1)+subFig(j).YLim(2); % Finds total y range from axis
        Cdata_val= findobj(subFig(j),'-property','CData'); % Loads color data value, creates x,y matrix (range depends on image steps)
        xi=round(xc*length(Cdata_val.CData)/xAxes); % adjusts input X value for range of available points in image
        yi=round(yc*length(Cdata_val.CData)/yAxes); % adjusts input Y value for range of available points in image
        
        xran=round(length(Cdata_val.CData(1,:))/2); % finds central X position for image
        yran=round(length(Cdata_val.CData(:,1))/2); % finds central Y position for image
        
        xv=xran+xi; % adds normalized user input to central X value     
        yv=yran+yi; % adds normalized user input to central Y value   
        
        spectra(i,j)=Cdata_val.CData(yv,xv); % finds spectral value for input 
    end
end
figure % plots and saves results 
subplot (2,3,1)
    plot(evRange,spectra(:,6))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(1));
subplot (2,3,2)
    plot(evRange,spectra(:,5))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(2));
subplot (2,3,3)
    plot(evRange,spectra(:,4))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(3));
subplot (2,3,4)
    plot(evRange,spectra(:,3))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(4));
subplot (2,3,5)
    plot(evRange,spectra(:,2))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(5));
subplot (2,3,6)
    plot(evRange,spectra(:,1))
    xlabel( 'eV (nm)' ); 
    ylabel( 'intensity (Au)' );
    title(name(6));
saveas(gcf,'extSpectra.fig')
saveas(gcf,'extSpectra.png')