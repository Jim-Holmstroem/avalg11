function [] = test_generator(filename,num_points) 
%
% Lets you generate a testfile from mouse input
% 
% saves it to filename
% num_points specifies the number of cities you wan't to plopout
%
%

points = zeros(2,num_points);

figure();
    axis manual;
    hold on;
    axis([0 15 0 10]);
    axis equal;
    
    for it = 1:num_points
        [x,y] = ginput(1);
        points(:,it)= [x;y];
        % plot(x,y,'x');
    end

    plot(points(1,:),points(2,:),'x');

    ginput(1);
    close all;


fid = fopen(filename,'w');

fprintf(fid,'%d\n',num_points);
for it = 1:num_points
    fprintf(fid,'%f %f\n',points(1,it),points(2,it));
end

st = fclose(filename);

