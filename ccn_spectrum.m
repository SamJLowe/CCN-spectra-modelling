function [CCN_spectrum] = ccn_spectrum(aerosol, CCN_spectrum)
%close all; clear all
%/////////////////////////////////////////////////////////////////////////%
%Function determines CCN number concentration (CCNC) as a function of 
%ambient supersaturation (samb) according to the critical supersaturations 
%stored in the aerosol sturcture (aerosol.crit_susat).
%---
%1. For given samb, identify dry aerosol size indices (k, k+1) which
%correspond to critical supersaturation (sc) that pinch samb.
%---
%2. (linearly, default)Interpolate to a dry radius between k and k+1, ract,
%corrsponding to sc = samb.
%---
%3. Use ract as part of the argument in complementary error function efrc
%to integrate the normally distributed log(r) distribution of aerosol conc.
%See e.g. Lowe et al. Nat. Comms vol. 10, Article number: 5214 (2019)
%methods. 
%---
% Note: Step 3 negates the need for explicit calculation of
% aerosol size distribution and any concerns about bin number (thus making 
% concentration interpolation in Lowe et al. (2016) unnecessary) or
% distribution limits in the integration step.
%---
%Indices:
%i - over supersaturation grid of CCN spectrum
%j - over aerosol modes
%k - over aerosol size bins
%/////////////////////////////////////////////////////////////////////////%

sc   = aerosol.crit_susat; rdry = aerosol.dry_radii;
ract = zeros(CCN_spectrum.nbins, aerosol.nmode);
%[naero] = aerosol_size_dist(aerosol);                                     %Can be used to sanity check erfc method against basic sum/integration

for i = 1:CCN_spectrum.nbins
    samb = CCN_spectrum.supersaturation(i);
    for j = 1:aerosol.nmode
        for k = 1:aerosol.nbins-1
            if samb < sc(j,k) && samb > sc(j,k+1)
                ract(i,j) = interp1([sc(j,k+1) sc(j,k)],...
                                    [rdry(k+1) rdry(k)],samb);                
                %CCNC_2(i,j) = sum(naero(j,rdry(2:end)>ract(i,j)))         %Can be used to sanity check erfc method against basic sum/integration                                      
            else 
            end
        end
    end
end

%use complementary error function erfc to integrate unitary normal distribution in
%log(r)-log(R)/log(GSD) from ract -> infinity. 
ract = ract';
x = (log(aerosol.modal_pars(:,2)*1e-9) - log(ract))  / sqrt(2) ./ log(aerosol.modal_pars(:,3)); 
CCNC = 0.5 * aerosol.modal_pars(:,1) .* (2 - erfc(x));

CCN_spectrum.ract = ract;
CCN_spectrum.CCNC = CCNC; 
%CCN_spectrum.CCNC = CCNC_2'; 

