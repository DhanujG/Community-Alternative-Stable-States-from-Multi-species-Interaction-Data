%%recreate Fig2e from Gonze et al. (2017).

%%set number of species;
S = 3;

%%set initial community abundances. These definitely matter and whoever
%%dominates depends on these. See Fig2 for initial densities for each
%%community state
x = zeros(S,1);
x(1) = 0.3;
x(2) = 0.1;
x(3) = 0.2;

%%construct K matrix for inhibition 
K = zeros(S);
inhibition_value = 0.1;

for ii = 1:S
    for jj = 1:S
        if ii == jj
            K(ii,jj) = 0;
        else
        K(ii,jj) = inhibition_value;
        end
    end
end   

%%define growth and death rates for simplest case (in future when
%%expanding for more than 3 species, we would write a more systematic way of defining these)
b = zeros(S,1);
k = zeros(S,1);
b(1) = 1;
b(2) = 0.95;
b(3) = 1.05;
k(1) = 1;
k(2) = 1;
k(3) = 1;

%%set number of timesteps for first iteration
timesteps = 100;

figure(1)
subplot(1,5,1)
[after_establishment] = Gonze_2017_calc_dynamics(x,b,k,K,timesteps);
title('Establishment of Community (b_1 = 1)')
ylabel('X_1,X_2,X_3')
species_names(1,1) = {'X_1'};
species_names(2,1) = {'X_2'};
species_names(3,1) = {'X_3'};
legend(species_names)

%%set perturbation (changing b_1 from Gonze et al.) and number of timesteps for second iteration
b(1) = 0.2;
timesteps = 30;

subplot(1,5,2)
[perturbation] = Gonze_2017_calc_dynamics(after_establishment,b,k,K,timesteps);
title('Perturbation (b_1 = 0.2)')

%%set back to initial and number of timesteps for second iteration
b(1) = 1;
timesteps = 100;

subplot(1,5,3)
[return_to_initial] = Gonze_2017_calc_dynamics(perturbation,b,k,K,timesteps);
title('Return from perturbation (b_1 = 1)')

%%set to stronger perturbation and number of timesteps for second iteration
b(1) = 4.5;
timesteps = 120;

subplot(1,5,4)
[stronger_perturbation] = Gonze_2017_calc_dynamics(return_to_initial,b,k,K,timesteps);
title('Stronger perturbation (b_1 = 4.5)')

%%set to initial and number of timesteps for second iteration
b(1) = 1;
timesteps = 200;

subplot(1,5,5)
[after_stronger] = Gonze_2017_calc_dynamics(stronger_perturbation,b,k,K,timesteps);
title('Return to initial (b_1 = 1)')