%% Refractive index data load for GUI ptot
function [n_loc,k_loc, enei_loc,n_loc2,k_loc2, enei_loc2, n_loc3, k_loc3, enei_loc3] = refPlot()
    % Load input variables  
    mat = dir('tmp/*.mat');
    for q = 1:length(mat) 
        load(sprintf(['tmp/' mat(q).name])); 
    end 
    nktmp=load(sprintf(['MNPBEM17/material/@epstable/' eps_part]));
    nktmp2=load(sprintf(['MNPBEM17/material/@epstable/' eps_part2]));
    nktmp3=load(sprintf(['MNPBEM17/material/@epstable/' eps_part3]));
    %% Refractive index Init
    n_loc=nktmp(:,2);
    k_loc=nktmp(:,3);
    enei_loc=nktmp(:,1);
    n_loc2=nktmp2(:,2);
    k_loc2=nktmp2(:,3);
    enei_loc2=nktmp2(:,1);
    n_loc3=nktmp3(:,2);
    k_loc3=nktmp3(:,3);
    enei_loc3=nktmp3(:,1);
end