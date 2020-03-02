% Put in folder with CL maps to extract from
% Verify only map .fig files are in the folder, after running remove produced
% .fig image to rerun
clear all
clc
xc=0; %input X value nm
yc=160; %input Y value nm

name=["CL", "CL p", "CL s", "CL_a", "CL_a p", "CL_a s"];
imagefiles = dir(sprintf('*.fig'));  % Find files   

for i=1:length(imagefiles)
    evRangeStr=imagefiles(i).name(8:end-6);
    evRange(i)=str2num(evRangeStr);
    fig=openfig(imagefiles(i).name, 'invisible'); % Loads new graphic
    subFig = findobj( gcf, 'Type', 'Axes' );
    for j=1:6
        xAxes=-subFig(j).XLim(1)+subFig(j).XLim(2);
        yAxes=-subFig(j).YLim(1)+subFig(j).YLim(2);
        Cdata_val= findobj(subFig(j),'-property','CData');
        xi=round(xc*length(Cdata_val.CData)/xAxes);
        yi=round(yc*length(Cdata_val.CData)/yAxes);
        
        xran=round(length(Cdata_val.CData(1,:))/2);
        yran=round(length(Cdata_val.CData(:,1))/2);
        
        if xi>2   
            xv=xran+xi;
        else
            xv=xran-xi;
        end
        if yi>2   
            yv=yran+yi;
        else
            yv=yran-yi;
        end        
        spectra(i,j)=Cdata_val.CData(yv,xv);
    end
end
figure
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