function hss_RunAllTests
%RUNALLTESTS Run all the unit tests.

addpath ../

hss_TestCreation;

hss_TestOperations;

hss_TestLyapunov;

rmpath ../

end
