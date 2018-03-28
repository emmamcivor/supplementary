function S_int=volume_integral_S(vol_int_S_z,vol_int_S_rt,Sub_conc_IC,particular_soln_S_all,r_s,theta,z_s)

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


% zero_mu=1:mu_max_s;
% mu=(2*zero_mu-1)*pi*0.5/(L2-L1);
% ef_z_m=cos((L2-z_s')*mu);
% ef_z_prime_int=sin(mu'*(L2-L1))./mu';
% 
% ef_z_z_prime_m=ef_z_m*diag(exp(-mu.^2*De*dt))*ef_z_prime_int;
% 
% Int_z_GF=2*ef_z_z_prime_m/(L2-L1);
% 
%             s0_z_min=min(min(S));
%             s0_z_max=max(max(S));
%             d=s0_z_max-s0_z_min;
%             edges=s0_z_min:d:s0_z_max;
%             if d==0
%                 S_int=S;
%             else
%                 [N_sub,edges_sub]=histcounts(S,edges);
% 
%                 background_conc_height_sub=max(N_sub);
%                 backround_conc_bin_sub=N_sub>=background_conc_height_sub;
%                 background_conc_sub=min(edges_sub(backround_conc_bin_sub));
%                 S_int=vol_int_S_z*(S-background_conc_sub)*vol_int_S_rt+Int_z_GF*background_conc_sub*ones(1,length(r_s)*length(theta))+particular_soln_S_all;
%                 
%             end

        
   

end