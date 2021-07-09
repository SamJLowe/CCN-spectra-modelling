%/////////////////////////////////////////////////////////////////////////%
%Top-level script for CCN spectra modelling using a suite of Köhler 
%formulations. Stockholm University (May, 2020). lowejsam@gmail.com for
%queries. 
%/////////////////////////////////////////////////////////////////////////%
close all; 
clear all; 
clc
tic

%supersaturations to calculate CCN concentration on
CCN_spectrum.nbins           = 100;
CCN_spectrum.supersaturation = linspace(0.01, 1.5, CCN_spectrum.nbins);

%setup aerosol data structure----------------------------------------------
aerosol.nmode      = 2;      %number of modes. each mode my have a different chemical makeup
aerosol.nbins      = 300;    %number of bins for size distribution
aerosol.rdry_min   = 10;      %minimum dry radius [nm]
aerosol.rdry_max   = 500;    %maximum dry radius [nm]

%Kohler model suite--------------------------------------------------------
%model = 'approx classic';   %[1] See e.g. pruppacher and klett 
model = 'classic';          %[2] See e.g. pruppacher and klett 
%model = 'Ovadnevaite';      %[3] See Ovadneviate et al. Nature 546, 637–641 supplementary info (2017).
%model = 'Ruehl';            %[4] See Ruehl et al. Science Vol. 351, Issue 6280, pp. 1447-1450 (2016).
%model = 'kappa';            %[5] See Petters and Kreidenweis Atmos. Chem. Phys., 7, 1961–1971, (2007).
%--------------------------------------------------------------------------

%set aerosol size distribution and composition
[aerosol]      = configure_aerosol(aerosol);
%calculate critical point of all aerosols according to chosen model
[aerosol]      = kohler_suite(model, aerosol);

% plot critical point------------------------------------------------------
figure(1); hold on
for i = 1:aerosol.nmode
plot(aerosol.crit_susat(i,:),aerosol.dry_radii*1e9,'.')
end
set(gca,'Yscale','log','Xscale','log'); xlabel('sc'); ylabel('rc [nm]')

%Determine CCN spectrum----------------------------------------------------
[CCN_spectrum] = ccn_spectrum(aerosol, CCN_spectrum);
toc
%Plot CCN spectrum and ract------------------------------------------------
figure(2)
    subplot(1,2,1)
        plot(CCN_spectrum.supersaturation, CCN_spectrum.ract(1,:)*1e9,'b.'); hold on
        plot(CCN_spectrum.supersaturation, CCN_spectrum.ract(2,:)*1e9,'r.')
        xlabel('susat [%]');ylabel('r_{act} [nm]');

    subplot(1,2,2)
        plot(CCN_spectrum.supersaturation, CCN_spectrum.CCNC(1,:),'b'); hold on
        plot(CCN_spectrum.supersaturation, CCN_spectrum.CCNC(2,:),'r')
        plot(CCN_spectrum.supersaturation, CCN_spectrum.CCNC(1,:) + CCN_spectrum.CCNC(2,:),'k')
        xlabel('susat [%]');ylabel('CCN conc [# cm^{-3}]');

%save output---------------------------------------------------------------
save('CCN_spectrum.mat','CCN_spectrum')

%Fin-----------------------------------------------------------------------