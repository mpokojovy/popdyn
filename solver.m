function [pop_m, pop_f] = solver(A_m, A_f, N_m, N_f, t_lat, eval_pts_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       US population dynamics simulation between 2001 and 2011           %
%      (C) Michael Pokojovy                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global maternity_data;
    global data_2001;
    global weighted_data_2001;
    global survival_data;
    global immigration_data;
    global weighted_immigration_data;
    
    load_population_statistics;
    
    prepare_maternity_data;
    prepare_init_data;
    prepare_immigration_data;

    h_m = A_m/(N_m - 1);
    h_f = A_f/(N_f - 1);
    
    AL_m = linspace(0, A_m, N_m);
    AL_f = linspace(0, A_f, N_f);
    ALI_m = AL_m(2:end);
    ALI_f = AL_f(2:end);
    
    e_m = ones(N_m - 1, 1)/h_m;
    e_f = ones(N_f - 1, 1)/h_f;
    
    A_mm = spdiags([e_m, -e_m], -1:0, N_m - 1, N_m - 1);
    A_ff = spdiags([e_f, -e_f], -1:0, N_f - 1, N_f - 1);
    A_mf = sparse(N_m - 1, N_f - 1);
    A_fm = sparse(N_f - 1, N_m - 1);
    A_mf(1, :) = A_mf(1, :) + maternity(0, ALI_f);
    A_ff(1, :) = A_ff(1, :) + maternity(1, ALI_f);
    
    A = [A_mm A_mf;
         A_fm A_ff];
     
    B = [zeros(1, N_f - 1) h_f*maternity(0, ALI_m)
         zeros(1, N_m - 1) h_f*maternity(1, ALI_f)];
     
    u0_m = init_data(0, ALI_m);
    u0_f = init_data(1, ALI_f);
    
    u0 = [u0_m; u0_f];
    
    sI = speye(N_m - 1 + N_f - 1);
    
    pop_m = zeros(N_m, length(eval_pts_ind));
    pop_f = zeros(N_f, length(eval_pts_ind));
    
    if length(find(eval_pts_ind == 1)) > 0
        pop_m(:, 1) = init_data(0, AL_m).*surv_prob(0, AL_m);
        pop_f(:, 1) = init_data(1, AL_f).*surv_prob(1, AL_f);
    end    
    
    for i = 1:length(t_lat) - 1
        dt = t_lat(i + 1) - t_lat(i);
        
        u0 = (sI/dt - 0.5*A)\((sI/dt + 0.5*A)*u0 + [imm_data(0, ALI_m); imm_data(1, ALI_f)]);
        
        ind = find(eval_pts_ind == i+1);
        
        if length(ind) > 0
            b = B*u0;
            
            u0_m = u0(1:N_m - 1);
            u0_f = u0(N_m:end);
            
            pop_m(:, ind) = [b(1); u0_m];
            pop_f(:, ind) = [b(2); u0_f];
            
            pop_m(:, ind) = pop_m(:, ind).*surv_prob(0, AL_m);
            pop_f(:, ind) = pop_f(:, ind).*surv_prob(1, AL_f);
        end
    end
    
    function m = maternity(sex, age)
        % male: sex == 0, female: sex == 1
        m = interp1(maternity_data(:, 1), maternity_data(:, 2 + sex), age);
    end

    function u = init_data(sex, age)
        u = interp1(weighted_data_2001(:, 1), weighted_data_2001(:, 2 + sex), age)';
    end

    function p = surv_prob(sex, age)
        p = interp1(survival_data(:, 1), survival_data(:, 2 + sex), age)';
    end

    function i = imm_data(sex, age)
        i = interp1(weighted_immigration_data(:, 1), weighted_immigration_data(:, 2 + sex), age)';
    end
    
    function prepare_maternity_data
        global birth_data;  
        
        maternity_data = birth_data;
        maternity_data(:, 2:3) = birth_data(:, 2:3).*survival_data(:, 2:3);
    end

    function prepare_init_data
        weighted_data_2001 = data_2001;
        weighted_data_2001(:, 2:3) = data_2001(:, 2:3)./survival_data(:, 2:3);
    end

    function prepare_immigration_data
        weighted_immigration_data = immigration_data;
        weighted_immigration_data(:, 2:3) = immigration_data(:, 2:3)./survival_data(:, 2:3);
    end
end