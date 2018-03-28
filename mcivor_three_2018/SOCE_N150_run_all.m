%% file with parameters to change in clustered configuration

%Parameters to change:
    distance_orai_SERCA=30;
    SERCA_choice=2;

% Simulate calcium concentrations    
    r_CRAC=30;
    NT=150;
    r_SERCA=r_CRAC+distance_orai_SERCA;
    length_t=1000;
    
	orai=[r_CRAC, r_SERCA]

    %% SERCA pump parameters (Lytton 1992)
    % SERCA2b
    Kmf_SERCA2b=0.27;
    n_H_SERCA2b=1.7;
    Vmax_ions_SERCA2b=36; %Ca2+ per second (according to Hogan 2015) 
    
    %SERCA2a
    Kmf_SERCA2a=0.38;
    n_H_SERCA2a=2.2;
    Vmax_ions_SERCA2a=72;  
    
    %% choose which SERCA pump is used
    if SERCA_choice==1
        disp('SERCA2a')
        Kmf=Kmf_SERCA2a; 
        n_H=n_H_SERCA2a;    
        n_avogadro=6.022e23;
        Vmax=(Vmax_ions_SERCA2a/n_avogadro)*1e6; % (36 ca ions/ avogradros number (ions per mol) * 1e6 = X micro moles per second)
    elseif SERCA_choice==2
        disp('SERCA2b')
        Kmf=Kmf_SERCA2b; 
        n_H=n_H_SERCA2b;    
        n_avogadro=6.022e23;
        Vmax=(Vmax_ions_SERCA2b/n_avogadro)*1e6; % (36 ca ions/ avogradros number (ions per mol) * 1e6 = X micro moles per second)
    end
        
        %% preload GF (check they are in memory or load if not)
    run SOCE_J_S_Orai_SERCA_Gauss_dr0p2_preload_GF_integrators_NT150.m
    
    %% run simulations with chosen Orai-SERCA distribution and SERCA pump isoform (2a or 2b)
    run SOCE_J_S_Orai_SERCA_Gauss_dr0p2_preloaded_GFI_NT150.m
    
    
toc
