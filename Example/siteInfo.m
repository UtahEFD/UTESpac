%For use with UTESpac
% site specific script to load site information
% Bonneville Salt Flats Experiment

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in ascending order
info.sonicOrientation = [25];

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1];  % ES5

% enter orientation of tower relative to sonic head
info.tower = 150;

% tower elevation
info.siteElevation = 1285.6; % (m)

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatnated with getData.m
info.tableNames = {'EC150_2'};

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [20];  %[Hz]

% enter number of columns in each .csv table.
info.tableNumberOfColumns = [16];
