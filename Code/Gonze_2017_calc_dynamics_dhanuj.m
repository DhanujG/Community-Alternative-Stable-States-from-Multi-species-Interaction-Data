function [final_densities,timeseries] = Gonze_2017_calc_dynamics_dhanuj(x,b,k,K,num_time_steps)
%Gonze_2017_calc_dynamics
%   Author: Dhanuj Gandikota, James Tan    Last Modified: April 15, 2020
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



    function [SSS] = base(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 1 = N1
        %N_1 = x(1);
        %R_N1 = b(1);
        Comp_N12 = 0.1;
        K1 = 0.8;
        c_R1 = 0.5
        Alpha_R1M = 1.1;
        h_R = 0.4;


        %       Resource Species 2 = N2
        %N_2 = x(2);
        %R_N2 = b(2);
        Comp_N21 = 0.2;
        K2 = 0.7;
        c_R2 = 0.5
        Alpha_R2M = 1.1;

        %   	Mutualistic Species M
%         R_M = 
         c_M = 0.8;
         KM = 0.7;
         Alpha_MR = 1.1;
%         Alpha_MS;
%         h_S;
         h_M = 0.4;
%         M;
%         d_M;
         q_M = 0.5;
         B_s = 0.6;
         e_s = 0.2;




        %   	Consumer Species 1 = C1
        %C_1 = x(3);
        b_1 = 1;
        %a_1 = b(3);
        %Comp_C12;
        h_C1 = 1;
        K_C1 = 0.8;
        X_1 = 0.05;

        %       Consumer Species 2 = C2
%         C_2;
%         b_2;
%         a_2;
%         Comp_C21;
%         h_C2;
%         K_C2;
%         X_2;

        %%Species Growth Functions

        %   RESOURCE Species Growth Function


        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 )) + c_M * ((Alpha_SM * N_1)/(h_M + M)) - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) - ( (a_2*C_2)/(1 + (a_2*h_C2*N_1)) ) );
        %dN2 = N_2*( (R_N2) * (1 - ( (N_2 + Comp_N21 * N_1)/K2 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_2)) ) - ( (a_2*C_2)/(1 + (a_2*h_C2*N_2)) ) );

        %               without mutualism and only 1 consumer
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );
        %dN2 = N_2*( (R_N2) * (1 - ( (N_2 + Comp_N21 * N_1)/K2 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_2)) ) );
        
        %               without mutualism and two consumers
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) - ( (a_2*C_2)/(1 + (a_2*h_C2*N_1)) ));
        
         %               without mutualism and one consumers and one
         %               resource species
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );

        %   CONSUMER Species Growth Function

        %dC1 = C_1*( (b_1 * (a_1 * N_1 + a_1 * N_2))/(1 + (a_1*h_C1*N_1) + (a_1*h_C1*N_2)) * (1 - ( (C_1 + Comp_C12 * C_2)/K_C1 )) + c_M * ((Alpha_SM * C_1)/(h_M + M)) - X_1);
        %dC2 = C_2*( (b_2 * (a_2 * N_1 + a_2 * N_2))/(1 + (a_2*h_C2*N_1) + (a_2*h_C2*N_2)) * (1 - ( (C_2 + Comp_C21 * C_1)/K_C2 )) - X_2);

        %               without mutualism one consumer
        %dC1 = C_1*( (b_1 * (a_1 * N_1 + a_1 * N_2))/(1 + (a_1*h_C1*N_1) + (a_1*h_C1*N_2)) * (1 - ( (C_1)/K_C1 ))  - X_1);
        
        %               without mutualism two consumer one resource
        %dC1 = C_1*( (b_1 * (a_1 * N_1 ))/(1 + (a_1*h_C1*N_1) ) * (1 - ( (C_1)/K_C1 ))  - X_1);
        %dC2 = C_2*( (b_2 * (a_2 * N_1 ))/(1 + (a_2*h_C2*N_1)) * (1 - ( (C_2)/K_C2 ))  - X_2);
        
        %               middle consumer and super consumer
        %dC1 = C_1*( (b_1 * (a_1 * N_1))/(1 + (a_1*h_C1*N_1)) * (1 - ( (C_1)/K_C1 )) - ( (a_2*C_2)/(1 + (a_2*h_C2*C_1)) )  - X_1);
        %dC2 = C_2*( (b_2 * (a_2 * C_1 ))/(1 + (a_2*h_C2*C_1)) * (1 - ( (C_2)/K_C2 ))  - X_2);

        %   Mutualistic Species Growth Function

        %dM = M * (R_M + c_M *((Alpha_MS* N_1)/(h_S + N_1)) - q_M * ((B_s * N_1)/(e_s + M)) - (d_M * M));
        %dM = M * (R_M + c_M *((Alpha_MS* C_1)/(h_S + C_1)) - q_M * ((B_s * C_1)/(e_s + M)) - (d_M * M));
        
        
        MRR = zeros(S,1);
        
        %      Consumer, Resource, Resource
        MRR(2) = x(2)*( (b(2)) * (1 - ( (x(2) + Comp_N12 * x(3))/K1 ))   + c_R1 * ((Alpha_R1M * x(2))/(h_M + x(1))));
        MRR(3) = x(3)*( (b(3)) * (1 - ( (x(3) + Comp_N21 * x(2))/K2 ))   + c_R2 * ((Alpha_R2M * x(3))/(h_M + x(1))));
        MRR(3) = x(1) * (((b(1)) * (1 - ( (x(1))/KM ))) + c_M *((Alpha_MR* x(2))/(h_R1 + x(2))) + c_M *((Alpha_MR* x(3))/(h_R2 + x(3))) - q_M * ((B_s * x(2))/(e_s + x(1))) - q_M * ((B_s * x(3))/(e_s + x(1))));
     end


    function [SRC] = map2a(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 1 = N1
        %N_1 = x(1);
        %R_N1 = b(1);
        Comp_N12 = 0.1;
        K1 = 0.8;


        %       Resource Species 2 = N2
        %N_2 = x(2);
        %R_N2 = b(2);
        Comp_N21 = 0.2;
        K2 = 0.7;




        %   	Consumer Species 1 = C1
        %C_1 = x(3);
        b_1 = 1;
        %a_1 = b(3);
        h_C1 = 1;
        K_C1 = 0.8;
        X_1 = 0.05;



        %%Species Growth Functions

        %   Resource Species Growth Function



        %               without mutualism and only 1 consumer for N1
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );
        %dN2 = N_2*( (R_N2) * (1 - ( (N_2 + Comp_N21 * N_1)/K2 ))  );

        %   Consumer Species Growth Function

        %               without mutualism, one consumer, one resource
        %               consumed
        %dC1 = C_1*( (b_1 * (a_1 * N_1))/(1 + (a_1*h_C1*N_1)) * (1 - ( (C_1)/K_C1 ))  - X_1);

        
        S
        SRC = zeros(S,1);
        
        %      Consumer, Resource, Resource
        SRC(1) = x(1)*( (b(1)) * (1 - ( (x(1) + Comp_N12 * x(2))/K1 ))  - ( (b(3)*x(3))/(1 + (b(3)*h_C1*x(1))) ) );
        SRC(2) = x(2)*( (b(2)) * (1 - ( (x(2) + Comp_N21 * x(1))/K2 )) );
        SRC(3) = x(3)*( (b_1 * (b(3) * x(1) + b(3) * x(2)))/(1 + (b(3)*h_C1*x(1)) + (b(3)*h_C1*x(2))) * (1 - ( (x(3))/K_C1 ))  - X_1);
     end


    function [RRC] = map2b(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 1 = N1
        %N_1 = x(1);
        %R_N1 = b(1);
        Comp_N12 = 0.1;
        K1 = 0.8;


        %       Resource Species 2 = N2
        %N_2 = x(2);
        %R_N2 = b(2);
        Comp_N21 = 0.2;
        K2 = 0.7;




        %   	Consumer Species 1 = C1
        %C_1 = x(3);
        b_1 = 1;
        %a_1 = b(3);
        h_C1 = 1;
        K_C1 = 0.8;
        X_1 = 0.05;



        %%Species Growth Functions

        %   Resource Species Growth Function



        %               without mutualism and only 1 consumer
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );
        %dN2 = N_2*( (R_N2) * (1 - ( (N_2 + Comp_N21 * N_1)/K2 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_2)) ) );

        %   Consumer Species Growth Function

        %               without mutualism one consumer
        %dC1 = C_1*( (b_1 * (a_1 * N_1 + a_1 * N_2))/(1 + (a_1*h_C1*N_1) + (a_1*h_C1*N_2)) * (1 - ( (C_1)/K_C1 ))  - X_1);

        
        
        RRC = zeros(S,1);
        
        %      Consumer, Resource, Resource
        RRC(1) = x(1)*( (b(1)) * (1 - ( (x(1) + Comp_N12 * x(2))/K1 ))  - ( (b(3)*x(3))/(1 + (b(3)*h_C1*x(1))) ) );
        RRC(2) = x(2)*( (b(2)) * (1 - ( (x(2) + Comp_N21 * x(1))/K2 ))  - ( ( b(3)*x(3))/(1 + (b(3)*h_C1*x(2))) ) );
        RRC(3) = x(3)*( (b_1 * (b(3) * x(1) + b(3) * x(2)))/(1 + (b(3)*h_C1*x(1)) + (b(3)*h_C1*x(2))) * (1 - ( (x(3))/K_C1 ))  - X_1);
    end

    function [CCR] = map2c(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 1 = N1
        %N_1 = x(1);
        %R_N1 = b(1);
        K1 = 0.8;




        %   	Consumer Species 1 = C1
        %C_1 = x(3);
        b_1 = 1;
        %a_1 = b(3);
        Comp_C12 = 0.5;
        h_C1 = 1;
        K_C1 = 0.8;
        X_1 = 0.1;

        %       Consumer Species 2 = C2
        %C_2 = x(2)
        b_2 = 0.5;
        %a_2 = b(2)
        Comp_C21 = 0.8;
        h_C2 = 1;
        K_C2 = 1;
        X_2 = 0.05;


        
        %               without mutualism and two consumers
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) - ( (a_2*C_2)/(1 + (a_2*h_C2*N_1)) ));

        
        %               without mutualism two consumer one resource
        %dC1 = C_1*( (b_1 * (a_1 * N_1 ))/(1 + (a_1*h_C1*N_1) ) * (1 - ( (C_1)/K_C1 ))  - X_1);
        %dC2 = C_2*( (b_2 * (a_2 * N_1 ))/(1 + (a_2*h_C2*N_1)) * (1 - ( (C_2)/K_C2 ))  - X_2);

        
        
        CCR = zeros(S,1);
        
        %      Consumer, Resource, Resource
        CCR(1) = x(1)*( (b(1)) * (1 - ( (x(1))/K1 ))  - ( (b(3)*x(3))/(1 + (b(3)*h_C1*x(1))) )  - ( (b(2)*x(2))/(1 + (b(2)*h_C2*x(1))) )) ;
        CCR(2) = x(2)*( (b_2 * (b(2) * x(1)))/(1 + (b(2)*h_C2*x(1))) * (1 - ( (x(2) + Comp_C21 * x(3))/K_C2 ))  - X_2);
        CCR(3) = x(3)*( (b_1 * (b(3) * x(1)))/(1 + (b(3)*h_C1*x(1))) * (1 - ( (x(3) + Comp_C12 * x(2))/K_C1 ))  - X_1);
     end
 
    function [CCRR] = map2d(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 1 = N1
        %N_1 = x(1);
        %R_N1 = b(1);
        %Comp_N12 = 0.1;
        K1 = 0.8;






        %   	Consumer Species 1 = C1
        %C_1 = x(3);
        b_1 = 1;
        %a_1 = b(3);
        h_C1 = 1;
        K_C1 = 0.8;
        X_1 = 0.05;

        %       Consumer Species 2 = C2
        %C_2 = x(2)
        b_2 = 0.5;
        %a_2 = b(2)
        h_C2 = 1;
        K_C2 = 1;
        X_2 = 0.05;

        %%Species Growth Functions

        %   RESOURCE Species Growth Function


        
         %               without mutualism and one consumers and one
         %               resource species
        %dN1 = N_1*( (R_N1) * (1 - ( (N_1)/K1 ))  - ( (a_1*C_1)/(1 + (a_1*h_C1*N_1)) ) );

        %   CONSUMER Species Growth Function

       
        
        %               middle consumer and super consumer
        %dC1 = C_1*( (b_1 * (a_1 * N_1))/(1 + (a_1*h_C1*N_1)) * (1 - ( (C_1)/K_C1 )) - ( (a_2*C_2)/(1 + (a_2*h_C2*C_1)) )  - X_1);
        %dC2 = C_2*( (b_2 * (a_2 * C_1 ))/(1 + (a_2*h_C2*C_1)) * (1 - ( (C_2)/K_C2 ))  - X_2);

 
        
        
        CCRR = zeros(S,1);
        
        %      Consumer, Resource, Resource
        CCRR(1) = x(1)*( (b(1)) * (1 - ( (x(1) )/K1 ))  - ( (b(2)*x(2))/(1 + (b(2)*h_C1*x(1))) ) );
        CCRR(2) = x(2)* (b_1 * (b(3) * x(1))/(1 + (b(3)*h_C1*x(1))) * (1 - ( (x(2))/K_C1 )) - ( (b(2)*x(3))/(1 + (b(2)*h_C2*x(2))) )   - X_1);
        CCRR(3) = x(3)*( (b_2 * (b(2) * x(2))/(1 + (b(2)*h_C2*x(2))) * (1 - ( (x(3))/K_C2 )) - X_2));
    end

    function [MRR] = map3a(~,x)

        %%Initialize these population functions

        %   Variables

        %       Resource Species 2 = N1
        %N_1 = x(2);
        %R_N1 = b(2);
        %Comp_N12 = 1.438187;
        Comp_N12 = 1.55619;
        K1 = 0.8;
        c_R1 = k(3);
        Alpha_R1M = 0.6;
        h_R = 0.3;


        %       Resource Species 3 = N2
        %N_2 = x(3);
        %R_N2 = b(3);
        Comp_N21 = 1;
        K2 = 0.8;
        c_R2 = k(4);
        Alpha_R2M = 0.6;

        %   	Mutualistic Species M
%         M = x(1)
  %       R_M = b(1)
         c_M1 = k(1);
         c_M2 = k(2);
         KM = 0.3;
         Alpha_MR = 0.6;
         %Alpha_MR = 0.6;
         h_M = 0.3;
         q_M1 = 0.5;
         q_M2 = 0.5;
         B_s = 0.2;
         e_s = 0.3;



        %%Species Growth Functions

        %   RESOURCE Species Growth Function


        %dN1 = N_1*( (R_N1) * (1 - ( (N_1 + Comp_N12 * N_2)/K1 )) + c_M * ((Alpha_SM * N_1)/(h_M + M)));


        %   Mutualistic Species Growth Function

        %dM = M * ((R_M * (1 - ( (M)/KM )))+ c_M *((Alpha_MS* N_1)/(h_S + N_1)) - q_M * ((B_s * N_1)/(e_s + M)) );
       
        
        MRR = zeros(S,1);
        
        %      Consumer, Resource, Resource
        MRR(2) = x(2)*( ((b(2)) *  ((K1 - x(2) - (Comp_N12 * x(3)))/K1 ))   + c_R1 * ((Alpha_R1M * x(1))/(h_M + x(1))));
        MRR(3) = x(3)*( ((b(3)) *  ((K2 - x(3) - (Comp_N21 * x(2)))/K2 ))   + c_R2 * ((Alpha_R2M * x(1))/(h_M + x(1))));
        MRR(1) = x(1) * (  ((b(1)) * (1 - ( (x(1))/KM ))) + c_M1 *((Alpha_MR* x(2))/(h_R + x(2))) + c_M2 *((Alpha_MR* x(3))/(h_R + x(3))) - q_M1 * ((B_s * x(2))/(e_s + x(1))) - q_M2 * ((B_s * x(3))/(e_s + x(1))));
     end 



%solve numerically
solution = ode15s(@map3a,[0 num_time_steps],x);

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