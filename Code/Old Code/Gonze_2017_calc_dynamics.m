function [final_densities,timeseries] = Gonze_2017_calc_dynamics(x,b,k,K,num_time_steps)
%Gonze_2017_calc_dynamics
%   Author: James Tan    Last Modified: April 12, 2020
%
%   [final_densities,timeseries] = Gonze_2017_calc_dynamics(x,b,k,K,num_tim
%   e_steps) solves the dynamical equations from
%   Gonze et al. 2017 toy model system and plots them on axes 
%
%   Use Gonze_2017_experiment.m to run the experiment to matching Gonze et al. 2017's 
%   Fig 2e, where the hysteresis dynamics are observed
%
%   Inputs: 
%   x is a vector of the initial community composition
%   b is a vector of the growth rate of each species
%   k is a vector of the death rate of each species
%   K is a matrix where K_ij is the inhibition of j on i, as defined by
%   Gonze et al.
%   num_time_steps is a salar that defined how long to run the simulation
%   for
%
%   Outputs: 
%   final_densities is a Sx1 vector with the final density of each arget
%       species at the end of the simulation.
%   timeseries is a matrix where each column of the first row is the 
%       timestep of the simulation and subsequent rows are the abundance of the
%       target species at each timestep. 

%determine species richness
S = numel(x);

%%set Hill coefficient
n = 2;

    % The function to numerically solve the set of coupled ordinary
    % differential equations defining the dynamics
    function [dxdt] = dynamics(~,x) 

        % calculate hill functions for each species
        
        hill = zeros(S);
        f = zeros(S,1);
        
        %%fill matrix with values of hill function values for each
        %%individual species where the value of the matrix at index i,j is
        %%the value of the Hill function of j on i. i on i has no
        %%inhibition
        for i = 1:S
            for j = 1:S
                if i == j
                    hill(i,j)=0;
                else
                    hill(i,j) = K(i,j)^n/(K(i,j)^n+x(j)^n);
                end
            end
        end
        
        %determine the f function using the values from above
        for i = 1:S
            f(i) = 1;
            for j = 1:S
               if i == j
                   f(i)=f(i);
               else
                   f(i) = f(i)*hill(i,j);
               end
            end
        end

        % calulate rate of change of densities; these can be vectorized 
        % for faster computation
        dxdt = zeros(S,1);
        
        for i = 1:S
            dxdt(i) = x(i)*(b(i)*f(i) - k(i)*x(i));
        end
    end


%     function [NNC] = wooster(~,x)
% 
% %       Resource Species 1 = N1
% 
%         Comp_N12 = 0.1;
%         K1 = 0.8;
% 
% 
% %       Resource Species 2 = N2
% 
%         Comp_N21 = 0.2;
%         K2 = 0.7;
% 
% 
% 
% 
% %   	Consumer Species 1 = C1
% 
%         b_1 = 1;
% 
%         h_C1 = 1;
%         K_C1 = 0.8;
%         X_1 = 0.05;
% 
% %       Consumer Species 2 = C2
%         %C_2;
%         %b_2;
%         %a_2;
%         %Comp_C21;
%         %h_C2;
%         %K_C2;
%         %X_2;
% 
% %               without mutualism and only 1 consumer
% %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );
% %dN2 = N_2*( (R_N2) * (1 - ( (N_2 + Comp_N21 * N_1)/K2 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_2)) ) );
% 
% 
% %               without mutualism one consumer
% %dC1 = C_1*( (b_1 * (a_1 * N_1 + a_1 * N_2))/(1 + (a_1*h_C1*N_1) + (a_1*h_C1*N_2)) * (1 - ( (C_1)/K_C1 ))  - X_1);
% 
% 
% %NNC(1) = x(1)*( (b(1)) * (1 - ( (x(1) + Comp_N12 * x(2))/K1 ))  - ( (0*b(3)*x(3))/(1 + (b(3)*h_C1*x(1))) ) );
% %NNC(2) = x(2)*( (b(2)) * (1 - ( (x(2) + Comp_N21 * x(1))/K2 ))  - ( (0* b(3)*x(3))/(1 + (b(3)*h_C1*x(2))) ) );
% %NNC(3) = x(3)*( (b_1 * (b(3) * x(1) + b(3) * x(2)))/(1 + (b(3)*h_C1*x(1)) + (b(3)*h_C1*x(2))) * (1 - ( (x(3))/K_C1 ))  - X_1);
% 
% 
% 
%         NNC = zeros(S,1);
%         NNC(1) = x(1)*( (b(1)) * (1 - ( (x(1) + Comp_N12 * x(2))/K1 ))  - ( (b(3)*x(3))/(1 + (b(3)*h_C1*x(1))) ) );
%         NNC(2) = x(2)*( (b(2)) * (1 - ( (x(2) + Comp_N21 * x(1))/K2 ))  - ( ( b(3)*x(3))/(1 + (b(3)*h_C1*x(2))) ) );
%         NNC(3) = x(3)*( (b_1 * (b(3) * x(1) + b(3) * x(2)))/(1 + (b(3)*h_C1*x(1)) + (b(3)*h_C1*x(2))) * (1 - ( (x(3))/K_C1 ))  - X_1);
%     end


%solve numerically
solution = ode15s(@dynamics,[0 num_time_steps],x);

% filter the solution so there's just one density for each species per 
% timestep
y = [x deval(solution,1:num_time_steps)];
timeseries = [(0:num_time_steps)' y'];

final_densities = y(:,end);

% create figure
% set species names for legend
% species_names(1,1) = {'X_1'};
% species_names(2,1) = {'X_2'};
% species_names(3,1) = {'X_3'};

% % choose species colors for consistency between figures
% colors = [236 128 46; 181 199 101; 124 136 63] ./256;
% colors = colors(species_names);

% plot the simulated population densities
h = plot(0:num_time_steps, y,'LineWidth',2); % DisplayName

% set plot labels
xlabel('Time')
% ylabel('X_1,X_2,X_3')
%legend(species_names)

% % specify display
% set(h,{'color'},num2cell(colors,2));
% set(gca,'yscale','log');
xlim([0 num_time_steps])
ylim([0 1])

end