%% Refractive index initialization
function epstab = refInit(Double_en,Cover_en, eps_env, eps_part, eps_str, eps_part2, eps_p, eps_list_en, eps_list_en2, substrate, substrate2, eps_p2, eps_part3, eps_list_en3, Source_op, qd_en)
if eps_list_en==0
    if eps_list_en2==0
        if Double_en==0 && substrate==0 %  Permittivity init
            if Cover_en==1 %  Permittivity init for cove layer
                epstab  = { epsconst( eps_env ), epsconst( eps_str ), epsconst( eps_p ) };
            else
                epstab = { epsconst( eps_env ), epsconst( eps_str) }; %  Table of dielectric function
            end
        else %  Permittivity init for double structures
            if substrate2==0
                epstab  = { epsconst( eps_env ), epsconst( eps_str ), epsconst( eps_p ) };
            else
                if eps_list_en3==0
                    epstab  = { epsconst( eps_env ), epsconst( eps_str ), epsconst( eps_p ), epsconst( eps_p2 ) };
                else
                    epstab  = { epsconst( eps_env ), epsconst( eps_str ), epsconst( eps_p ), epstable( eps_part3 ) };
                end
            end
        end
    else
        if Double_en==0 && substrate==0 %  Permittivity init
            if Cover_en==1 %  Permittivity init for cove layer
                epstab  = { epsconst( eps_env ), epsconst( eps_str ), epstable( eps_part2 ) };
            else
                epstab = { epsconst( eps_env ), epsconst( eps_str) }; %  Table of dielectric function
            end
        else %  Permittivity init for double structures
            if substrate2==0
                epstab  = { epsconst( eps_env ), epsconst( eps_str ), epstable( eps_part2 ) };
            else
                if eps_list_en3==0
                    epstab  = { epsconst( eps_env ), epsconst( eps_str ), epstable( eps_part2 ), epsconst( eps_p2 ) };
                else
                    epstab  = { epsconst( eps_env ), epsconst( eps_str ), epstable( eps_part2 ), epstable( eps_part3 ) };
                end
            end
        end
    end
else
    if eps_list_en2==0
        if Double_en==0 && substrate==0 %  Permittivity init
            if Cover_en==1 %  Permittivity init for cove layer
                epstab  = { epsconst( eps_env ), epstable( eps_part ), epsconst( eps_p ) };
            else
                epstab = { epsconst( eps_env ), epstable( eps_part ) }; %  Table of dielectric function
            end
        else %  Permittivity init for double structures
            if substrate2==0
                epstab  = { epsconst( eps_env ), epstable( eps_part ), epsconst( eps_p ) };
            else
                if eps_list_en3==0
                    epstab  = { epsconst( eps_env ), epstable( eps_part ), epsconst( eps_p ), epsconst( eps_p2 ) };
                else
                    epstab  = { epsconst( eps_env ), epstable( eps_part ), epsconst( eps_p ), epstable( eps_part3 ) };
                end
            end
        end
    else
        if Double_en==0 && substrate==0 %  Permittivity init
            if Cover_en==1 %  Permittivity init for cove layer
                epstab  = { epsconst( eps_env ), epstable( eps_part ), epstable( eps_part2 ) };
            else
                epstab = { epsconst( eps_env ), epstable( eps_part ) }; %  Table of dielectric function
            end
        else %  Permittivity init for double structures
            if substrate2==0
                epstab  = { epsconst( eps_env ), epstable( eps_part ), epstable( eps_part2 ) };
            else
                if eps_list_en3==0
                    epstab  = { epsconst( eps_env ), epstable( eps_part ), epstable( eps_part2 ), epsconst( eps_p2 ) };
                else
                    epstab  = { epsconst( eps_env ), epstable( eps_part ), epstable( eps_part2 ), epstable( eps_part3 ) };
                end
            end
        end
    end
end   
if qd_en==1 && Source_op==3
    if eps_list_en3==0
        epstab{end+1} = epsconst( eps_p2 );
    else
        epstab{end+1} = epstable( eps_part3 );
    end    
end
end