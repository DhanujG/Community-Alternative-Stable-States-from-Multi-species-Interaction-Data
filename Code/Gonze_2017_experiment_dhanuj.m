%%recreate Fig2e from Gonze et al. (2017).

%%set number of species;
S = 3;

%%set initial community abundances. These definitely matter and whoever
%%dominates depends on these. See Fig2 for initial densities for each
%%community state
x = zeros(S,1);
x(1) = 0.2;
x(2) = 0.1;
x(3) = 0.7;

%%construct K matrix for inhibition 
K = zeros(S);
% inhibition_value = 0.1;
% 
% for ii = 1:S
%     for jj = 1:S
%         if ii == jj
%             K(ii,jj) = 0;
%         else
%         K(ii,jj) = inhibition_value;
%         end
%     end
% end   

%%define growth and death rates for simplest case (in future when
%%expanding for more than 3 species, we would write a more systematic way of defining these)
b = zeros(S,1);
k = zeros(S + 1,1);
b(1) = 0.2;
b(2) = 0.2;
b(3) = 0.2;
k(1) = 0.3;
k(2) = 0.2;
k(3) = 0.5;
k(4) = 0.1;

%%set number of timesteps for first iteration
timesteps = 1000;

figure(1)
subplot(1,5,1)
[after_establishment] = Gonze_2017_calc_dynamics_dhanuj(x,b,k,K,timesteps);
title('Establishment of Community M - R - R ')
ylabel('M,R_1,R_2')
species_names(1,1) = {'Mutualism'};
species_names(2,1) = {'Resource_1'};
species_names(3,1) = {'Resource_2'};
legend(species_names)

%%set perturbation (changing b_1 from Gonze et al.) and number of timesteps for second iteration
k(1) = k(1)*1.1;
%b(3) = b(3)*2;
timesteps = 60;

subplot(1,5,2)
[perturbation] = Gonze_2017_calc_dynamics_dhanuj(after_establishment,b,k,K,timesteps);
title('Perturbation (c_M1 * 1.1)')

%%set back to initial and number of timesteps for second iteration
k(1) = k(1)/1.1;
%b(3) = b(3)/2;

subplot(1,5,3)
[return_to_initial] = Gonze_2017_calc_dynamics_dhanuj(perturbation,b,k,K,timesteps);
title('Return from perturbation ')

%%set to stronger perturbation and number of timesteps for second iteration
k(1) = k(1)/5;
%b(3) = b(3)*4;
timesteps = 240;

subplot(1,5,4)
[stronger_perturbation] = Gonze_2017_calc_dynamics_dhanuj(return_to_initial,b,k,K,timesteps);
title('Perturbation (c_M2 / 5)')

%%set to initial and number of timesteps for second iteration
k(1) = k(1)*5;
%b(3) = b(3)/2;
timesteps = 1000;

subplot(1,5,5)
[after_stronger] = Gonze_2017_calc_dynamics_dhanuj(stronger_perturbation,b,k,K,timesteps);
title('Return to initial')