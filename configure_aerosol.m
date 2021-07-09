function [aerosol] = configure_aerosol(aerosol)

%Aerosol size dsitribution description ------------------------------------
aerosol.dry_radii  = linspace(aerosol.rdry_min, aerosol.rdry_max,aerosol.nbins)*1e-9;      %bin radii [m]
aerosol.volumes    = 4 * pi * aerosol.dry_radii.^3 / 3;                    %aerosol volume [m3]
aerosol.modal_pars =  [265, 21.0, 1.45;...    %log-normal pars: number conc. [#/cm3], geometric mean radius [nm], geometric std
                       165, 82.5, 1.50];      %second mode...number of rows must = aerosol.nmode
              
%Chemical configuration---------------------------------------------------- 
%These variables are used in models [1-4] only-----------------------------
                      %model organic(succinic acid), black carbon, ammonium sulfate,
                      %sodium chlorid, ammonium nitrate
aerosol.comp_dens   = [1560, 2000, 1770, 2160, 1720;... %[kg/m3] elements = [compound, modes]
                       1560, 2000, 1770, 2160, 1720];
aerosol.molar_mass  = [118.09 , 12, 132, 58.44, 80.55;... %[kg/mol] elements = [compound, modes]
                       118.09 , 12, 132, 58.44, 80.55]*1e-3;
                   FORG = 0.2
aerosol.vol_fracs   = [FORG, 0.0, 1-FORG, 0.0, 0.0;...  %volume fraction of each compound                             
                       FORG, 0.0, 1-FORG, 0.0, 0.0];                                    
aerosol.comp_dissoc = [1, 0, 3, 2, 2;...               %ionic dissociation
                       1, 0, 3, 2, 2];

%These variables are used in model [3] only--------------------------------
%generic pars from ranges in ovadnevaite et al SI
aerosol.min_film_thickness  = [0.2e-9, 0.2e-9]; %[m] elements = modes
aerosol.surface_tension_org = [40e-3, 40e-3];   %[Nm-1] elements = modes
%--------------------------------------------------------------------------

%These variables are used in model [4] only--------------------------------
%WARNING: these parameters is compound specific, here is for succinic acid,
%more can be found in Ruehl et al., and organic properties above must be
%adjusted to the same compound
 aerosol.A0      = [76.9, 76.9]*1e-20; %[Area]= [m^2]   
 aerosol.C0      = [6e-7, 6e-7];  
 aerosol.m_sigma = [0.00104, 0.00104]*1e20; % [Jm-2m2]
 aerosol.surftens_min = [0.0354, 0.0354]; %[Jm-2]
%--------------------------------------------------------------------------

%These variables are used in model [5] only--------------------------------
%values based on nothing...
aerosol.kappa   = [.5, .5];
%--------------------------------------------------------------------------
end

