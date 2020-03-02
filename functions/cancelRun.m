%% Cancel function, Cancels running calculations
function cancelRun()        
    flag=load('tmp/cFlag.mat'); % Loads cancel flag 
    if flag.cFlag==1
        error('Program terminated by user')
    end
end