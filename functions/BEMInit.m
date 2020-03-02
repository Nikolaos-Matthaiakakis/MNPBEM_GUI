%% BEM initialization, sets solver options
function op = BEMInit(BEM_op, neg, intPoint, dIntPoint, Relcut)
    if BEM_op==0
        op = bemoptions( 'sim', 'ret', 'interp', 'curv', 'refine', intPoint, 'npol', dIntPoint, 'RelCutoff' , Relcut ); % Options for BEM (Retarded)
    elseif BEM_op==1
        op = bemoptions( 'sim', 'stat', 'interp', 'curv','nev', neg, 'refine', intPoint, 'npol', dIntPoint, 'RelCutoff' , Relcut ); % Options for BEM (Quasistatic)
    else % Options for BEM (Iterative)
        %  option structure
        op = bemoptions( 'sim', 'ret','interp', 'curv', 'refine', intPoint, 'npol', dIntPoint, 'RelCutoff' , Relcut  );
        %  add options for iterative solver
        op.iter = bemiter.options;
    end
end
