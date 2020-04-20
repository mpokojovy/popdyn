%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       US population dynamics simulation between 2001 and 2011           %
%      (C) Michael Pokojovy                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global data_2001 data_2011;

T = 10; % number of years from 2001

N_m = 100*12; % size of the age lattice for male population
N_f = 100*12; % size of the age lattice for female population

A_m = 100; % age lattice for male population
A_f = 100; % age lattice for female population

eval_pts_ind = 10*12;
t_lat = linspace(0, 10, eval_pts_ind); % time lattice

[pop_m, pop_f] = solver(A_m, A_f, N_m, N_f, t_lat, eval_pts_ind);

age = (0:100)';
pop_m = interp1(linspace(0, A_m, N_m), pop_m, age);
pop_f = interp1(linspace(0, A_f, N_f), pop_f, age);

%%%
figure(1);
plot_population([age pop_m pop_f]);
mtit('Simulated U.S. population structure in 2011', 'interpreter', 'latex', 'FontSize', 18);
%%
figure(2);
plot_population(data_2011);
mtit('Reported U.S. population structure in 2011', 'interpreter', 'latex', 'FontSize', 18);
%%

male_pop_size_2011_rep = sum(data_2011(:, 2));
female_pop_size_2011_rep = sum(data_2011(:, 3));
male_pop_size_2011_sim = sum(pop_m);
female_pop_size_2011_sim = sum(pop_f);

display(['U.S. male population in 2011']);
display(' ');
display(['Reported size ', num2str(round(male_pop_size_2011_rep))]);
display(['Simulated size ', num2str(round(male_pop_size_2011_sim))]);
display(['Relative discrepancy ', num2str(abs(male_pop_size_2011_sim - male_pop_size_2011_rep)/male_pop_size_2011_rep)]);

display(['Absolute L^1 error ', num2str(norm(pop_m - data_2011(:, 2), 1))]);
display(['Absolute L^2 error ', num2str(norm(pop_m - data_2011(:, 2), 2))]);
display(['Absolute L^inf error ', num2str(norm(pop_m - data_2011(:, 2), Inf))]);

display(['Relative L^1 error ', num2str(norm(pop_m - data_2011(:, 2), 1)/male_pop_size_2011_rep)]);
display(['Relative L^2 error ', num2str(norm(pop_m - data_2011(:, 2), 2)/male_pop_size_2011_rep)]);
display(['Relative L^inf error ', num2str(norm(pop_m - data_2011(:, 2), Inf)/male_pop_size_2011_rep)]);

display(' ');

display(['U.S. female population in 2011']);
display(' ');
display(['Reported size ', num2str(round(female_pop_size_2011_rep))]);
display(['Simulated size ', num2str(round(female_pop_size_2011_sim))]);
display(['Relative discrepancy ', num2str(abs(female_pop_size_2011_sim - female_pop_size_2011_rep)/female_pop_size_2011_rep)]);

display(['Absolute L^1 error ', num2str(norm(pop_f - data_2011(:, 3), 1))]);
display(['Absolute L^2 error ', num2str(norm(pop_f - data_2011(:, 3), 2))]);
display(['Absolute L^inf error ', num2str(norm(pop_f - data_2011(:, 3), Inf))]);

display(['Relative L^1 error ', num2str(norm(pop_f - data_2011(:, 3), 1)/female_pop_size_2011_rep)]);
display(['Relative L^2 error ', num2str(norm(pop_f - data_2011(:, 3), 2)/female_pop_size_2011_rep)]);
display(['Relative L^inf error ', num2str(norm(pop_f - data_2011(:, 3), Inf)/female_pop_size_2011_rep)]);

display(' ');