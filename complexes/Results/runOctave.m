% Test file
% disp('This is a test')

A = dlmread("bonds.txt");
B = dlmread("energy.txt");

C = A \ B;
dlmwrite("resCoef", C);

disp('Octave script finished')

