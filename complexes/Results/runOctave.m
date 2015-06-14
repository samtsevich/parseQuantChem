% Test file
% disp('This is a test')

A = dlmread("bonds.txt");
b = dlmread("energy.txt");

c = A \ b;
dlmwrite("resCoef", c);
dlmwrite("resEner", A*c);

disp('Octave script finished')

