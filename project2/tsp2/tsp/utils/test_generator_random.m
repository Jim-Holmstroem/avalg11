function [] = test_gernerator_random(filename,num_points)




points = 10*rand(2,num_points);

figure();
    axis manual;
    hold on;
    axis([0 15 0 10]);
    axis equal;

    plot(points(1,:),points(2,:),'x')
    ginput(1);
    close all;

fid = fopen(filename,'w');

fprintf(fid,'%d\n',num_points);
for it = 1:num_points
    fprintf(fid,'%f %f\n',points(1,it),points(2,it));
end

st = fclose(filename);



