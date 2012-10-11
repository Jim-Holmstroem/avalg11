function [] = plotter(cities_file,solution_file)
% Plots the cities and solution if it has the same format as in the problem: oldkattis:tsp
% 
% kattis input =   cities_file
% kattis output= solution_file
%
% cities_file format:
% <num of cities>
% <x_1> <y_1>
% <x_2> <y_2>
% ...
% <x_num> <y_num>
%
% solution_file format (all 0-indexed integers representing the city order):
% <s_1>
% <s_2>
% ...
% <s_num>
%
% If no argument is specified the default files ['cities.dat','solution.dat'] will be used. 
%


if(nargin<2)
    cities_file = 'cities.dat';
    solution_file = 'solution.dat';
end

fid = fopen(cities_file,'rt');
num = fscanf(fid,'%d\n',1);
[cities,count] = fscanf(fid,'%f %f\n',2*num);
fclose(fid);

fid = fopen(solution_file,'rt');
[solution,sol_count] = fscanf(fid,'%d',num);
fclose(fid);

cities = reshape(cities,2,num);
solution=solution+1; %make it 1-indexed
solution=[solution;solution(1)]; %the last distance is also a a part of the solution

figure();
	plot(cities(1,:),cities(2,:),'x');
	hold on;
	plot(cities(1,solution),cities(2,solution),'r')
	
end

