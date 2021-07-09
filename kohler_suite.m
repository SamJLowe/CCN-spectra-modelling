function aerosol = kohler_suite(model, aerosol)

%/////////////////////////////////////////////////////////////////////////%
%Function receives aerosol data structure and string for model name and
%returns critical point of köhler curve in the aerosol structure.
%/////////////////////////////////////////////////////////////////////////%

%some physical constants
rhow        = 1000;            %water density [kg/m3]
Mw          = 18e-3;           %water molar mass [kg/mol]
NA          = 6.0221409e23;    %Avagadros num.
T           = 285;             %temperature [K]
R           = 8.314;           %Universal gas constant [J K-1 mol-1]
surf_tens_w = 72.8e-3;         %surface tension of water droplet [Nm-1]
V_w         = Mw/rhow;

%for specifying wet raidii to loop over
wet_radius_end  = 10e-6;
wet_radius_bins = 300;

%Köhler A coefficent, used throughout
%A = 2 * Mw * surf_tens_w / R / T / rhow;
%--------------------------------------------------------------------------
switch model
    %----------------------------------------------------------------------
    case 'classic'
        Acoeff = 2 * Mw * surf_tens_w / R / T / rhow;
        for i = 1:aerosol.nmode
            factor(i) = sum(aerosol.comp_dens(i,:).*aerosol.comp_dissoc(i,:).*... %soluble Moles/volumes
                aerosol.vol_fracs(i,:)./aerosol.molar_mass(i,:));
            for j = 1:aerosol.nbins
                wet_radius = linspace(aerosol.dry_radii(j)*1.1,...
                    wet_radius_end, wet_radius_bins);
                for k = 1:wet_radius_bins
                    kelvin_term(k) = exp(Acoeff / wet_radius(k));
                    moles_water = 4 * pi / 3 * (wet_radius(k)^3 - aerosol.dry_radii(j)^3) / V_w;
                    soluble_moles = factor(i) * aerosol.volumes(j);                    
                    x(k) = moles_water / (moles_water + soluble_moles);
                end
                S = x.*kelvin_term;
                [aerosol.crit_susat(i,j), icrit] =(max(100 * (S-1)));
                aerosol.crit_radius(i,j) = wet_radius(icrit);
                %figure(1); plot(wet_radius,kelvin_term); hold on
            end
        end
        %------------------------------------------------------------------
    case 'approx classic'
        Acoeff =  3.3e-5 / T;
        Bcoeff   = zeros(aerosol.nmode,aerosol.nbins);
        for i = 1:aerosol.nmode
            factor(i) = sum(aerosol.comp_dens(i,:).*aerosol.comp_dissoc(i,:).*... %soluble Moles/volumes
                aerosol.vol_fracs(i,:)./aerosol.molar_mass(i,:));
            for j = 1:aerosol.nbins
                soluble_moles = aerosol.volumes(j) * factor(i);
                Bcoeff       = 4.3 * soluble_moles;
                aerosol.crit_radius(i,j) = (3 * Bcoeff / Acoeff)^0.5;
                aerosol.crit_susat(i,j) = 100 * (4 * Acoeff^3 / 27 / Bcoeff)^0.5;
            end
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
    case 'Ovadnevaite'
        aerosol.comp_dissoc(:,1) = 0; %insolbule organic
        surface_tension = zeros(aerosol.nmode,aerosol.nbins,wet_radius_bins);
        kelvin_term     = zeros(1,wet_radius_bins);
        x               = zeros(1,wet_radius_bins);
        factor          = zeros(1,aerosol.nmode);
        soluble_moles   = zeros(aerosol.nmode,aerosol.nbins);
        for i = 1:aerosol.nmode
            %model assumes all organic material is at surface in a film state
            film_volume = aerosol.volumes * aerosol.vol_fracs(i,1);
            %total soluble dry moles/volume in particle bulk phase
            factor(i) = sum(aerosol.comp_dens(i,:).*aerosol.comp_dissoc(i,:).*...
                aerosol.vol_fracs(i,:)./aerosol.molar_mass(i,:));
            for j = 1:aerosol.nbins
                %soluble moles of a particle in a given size bin
                soluble_moles(i,j) = (aerosol.volumes(j)-film_volume(j)) * factor(i);
                wet_radius = linspace(aerosol.dry_radii(j)*1.001,...
                    wet_radius_end, wet_radius_bins);
                for k = 1:wet_radius_bins
                    %minimum film volume for minimum film thickness
                    min_film_volume = 4 * pi / 3 *(wet_radius(k)^3 - ...
                        (wet_radius(k) - aerosol.min_film_thickness(i))^3);
                    %fractional surface coverage of organic film
                    surf_coverage   = min([film_volume(j)/min_film_volume, 1]);
                    %Surface coverage weighted surface tension
                    surface_tension(i,j,k) = surf_coverage * aerosol.surface_tension_org(i) +...
                        (1 - surf_coverage) * surf_tens_w;
                    %calculate kelvin term
                    kelvin_term(k) = exp(2 * Mw * surface_tension(i,j,k) / R / T / rhow /wet_radius(k));
                    %moles of water
                    nw = rhow / Mw * 4 * pi / 3 * wet_radius(k) ^ 3;
                    %mole fraction
                    x(k) = nw / (nw + soluble_moles(i,j));
                    
                end
                S = x .* kelvin_term;
                %                 figure(3); subplot(1,2,1); plot(wet_radius, squeeze(surface_tension(i,j,:))); hold on; set(gca,'Xscale','log')
                %                 figure(3); subplot(1,2,2); plot(wet_radius, S,'k'); hold on;
                %                 plot(wet_radius, x,'b')
                %                 plot(wet_radius, kelvin_term,'r')
                %                 set(gca,'Xscale','log')
                
                
                [aerosol.crit_susat(i,j), icrit] =(max(100 * (S-1)));
                aerosol.crit_radius(i,j) = wet_radius(icrit);
                aerosol.crit_surface_tension(i,j) = surface_tension(i,j,icrit);
            end
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    case 'Ruehl'
        for i = 1:aerosol.nmode
            V_org = aerosol.molar_mass(i,1)/aerosol.comp_dens(i,1); %*1e27
            factor(i) = sum(aerosol.comp_dens(i,:).*aerosol.comp_dissoc(i,:).*...
                aerosol.vol_fracs(i,:)./aerosol.molar_mass(i,:));
            for j = 1:aerosol.nbins
                wet_radius = linspace(aerosol.dry_radii(j)*1.001,...
                    wet_radius_end, wet_radius_bins);
                for k = 1:wet_radius_bins
                    
                    %solving the isotherm.--
                    seed_radius3 = 3 * aerosol.volumes(j) * (1-aerosol.vol_fracs(i,1)) / 4 / pi;          %associated with inorganic mass
                    A_iso = 3 *  V_org *wet_radius(k)^2 / ((aerosol.dry_radii(j)^3 - seed_radius3)*NA);   %actually A * fsurf
                    C_bulk = (aerosol.dry_radii(j)^3 - seed_radius3)*V_w/((wet_radius(k)^3)*V_org);       %actually C_bulk  / (1-fsurf)
                    %Compressed film isotherm to solve
                    iso = @(fsurf) (1-fsurf)*C_bulk - aerosol.C0(i) * ...
                        exp((aerosol.A0(i)^2-(A_iso/fsurf)^2)*aerosol.m_sigma(i)*NA/2/R/T);
                    %solver
                    f_surf(k) = fzero(iso,0.5);
                    
                    %----------------------
                    %bilayer
                    %A_min=aerosol.A0(i)+(aerosol.surftens_min(i)-0.0728)/aerosol.m_sigma(i);
                    %                      if A_iso/f_surf(k) < A_min
                    %                          f_surf(k) = pi * (wet_radius(k)*2)^2/(moles_org*A_min*NA);
                    %                      end
                    %----------------
                    
                    surface_tension(k) = max([surf_tens_w - (aerosol.A0(i) - A_iso/f_surf(k)) * aerosol.m_sigma(i), aerosol.surftens_min(i)]);
                    surface_tension(k) = min([surface_tension(k), surf_tens_w]);
                    moles_aero_total = factor(i)* aerosol.volumes(j);
                    moles_org = 4 * pi / 3 * (aerosol.dry_radii(j)^3 - seed_radius3) / V_org;             %total moles of organic
                    moles_inorg = moles_aero_total - moles_org;
                    moles_water = 4 * pi / 3 * (wet_radius(k)^3 - seed_radius3) / V_w;
                    moles_inorg = factor(i)*(aerosol.volumes(j));
                    
                    %moles fraction of water (water activity assuming ideality)
                    x(k) = moles_water / (moles_water + moles_org * (1-f_surf(k)) + moles_inorg);
                    kelvin_term(k) = exp(2 * Mw * surface_tension(k) / R / T / rhow /wet_radius(k));
                    
                end
                S = x .* kelvin_term;
                %                 figure(4)
                %                 subplot(4,1,1:2);plot(wet_radius*2,S,'.'); hold on; xlim([800e-9 2400e-9])
                %                 subplot(4,1,3);plot(wet_radius*2,f_surf,'.'); hold on; ylim([0 1]); xlim([800e-9 2400e-9])
                %                 subplot(4,1,4);plot(wet_radius*2,surface_tension,'.'); hold on; xlim([800e-9 2400e-9])
                
                
                
                [aerosol.crit_susat(i,j), icrit] =(max(100 * (S-1)));
                aerosol.crit_radius(i,j) = wet_radius(icrit);
                aerosol.crit_surface_tension(i,j) = surface_tension(icrit);
            end
        end
        %----------------------------------------------------------------------
        %----------------------------------------------------------------------
    case 'kappa'
        Acoeff = 2 * Mw * surf_tens_w / R / T / rhow;
        for i = 1:aerosol.nmode
            for j = 1:aerosol.nbins
                wet_radius = linspace(aerosol.dry_radii(j)*1.001,...
                    wet_radius_end, wet_radius_bins);
                for k = 1:wet_radius_bins
                    kelvin_term(k) = exp(Acoeff / wet_radius(k));
                    kappa_term(k) = (wet_radius(k)^3 - aerosol.dry_radii(j)^3) / ...
                        (wet_radius(k)^3 - aerosol.dry_radii(j)^3 * ...
                        (1 - aerosol.kappa(i)));
                end
                S = kelvin_term .* kappa_term;
                %figure(3); plot(wet_radius, S); hold on; set(gca,'Xscale','log')
                [aerosol.crit_susat(i,j), icrit] =(max(100 * (S-1)));
                aerosol.crit_radius(i,j) = wet_radius(icrit);
                
            end
        end
        %----------------------------------------------------------------------
    otherwise
        error('Incorrect model string')
end
