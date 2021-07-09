function [naero] = aerosol_size_dist(aerosol)

%aerosol = configure_aerosol();
rdry = aerosol.dry_radii;
N = aerosol.modal_pars(:,1); 
R = aerosol.modal_pars(:,2) * 1e-9; 
GSD = aerosol.modal_pars(:,3);

for i = 1:aerosol.nmode

dlogr = diff(log(aerosol.dry_radii));    
naero(i,:) = dlogr .* N(i) / sqrt(2*pi) / log(GSD(i)) .*...
exp(-1*(log(rdry(2:end)/R(i))).^2/ 2 / log(GSD(i))^2)
%plot(rdry(2:end), naero(i,:)); hold on; set(gca,'Xscale','log')

end


