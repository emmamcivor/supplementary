function S_int=volume_integral_S_flux_BC(vol_int_S_z,vol_int_S_rt,Sub_conc_IC,particular_soln_S_all,z_s,r_s,theta)

    S=Sub_conc_IC-particular_soln_S_all;
            
    background_conc_sub=zeros(1,length(z_s));
        for pp=1:length(z_s)
            s0_z_min=min(S(pp,:));
            s0_z_max=max(S(pp,:));
            d=(s0_z_max-s0_z_min)/100;
            [N_sub,edges_sub]=histcounts(S(pp,:));
            if d==0
                background_conc_sub(pp)=s0_z_min;
            else
                background_conc_height_sub=max(N_sub);
                backround_conc_bin_sub=N_sub>=background_conc_height_sub;
                background_conc_sub(pp)=min(edges_sub(backround_conc_bin_sub));
            end

        end

    background_conc_sub_remove=repmat(background_conc_sub,[length(r_s)*length(theta) 1])';
    S_int=vol_int_S_z*(S-background_conc_sub_remove)*vol_int_S_rt+particular_soln_S_all+vol_int_S_z*background_conc_sub_remove;

end