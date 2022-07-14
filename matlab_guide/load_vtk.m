% VTK_read - read from VTK file 
%
%   VTKData = VTK_read(filename,ext)
%
%   VTKData is a structure currently containing:
%       1. The vtk file headers
%       2. The grid type
%       3. Point data
%       4. Cell radius
%       5. Cell status
%       6. Cell concentration
%       7. fibres
%       8. phenot
%       9. Cell phenotype
%
%   Adapted from a script 
%   Copyright 2022 (c) Zhangxi Feng, University of New Hampshire

%% INSTRUCTIONS 
% This script is executed first before make_plots
% Copy the file "first_simulation_cells.1000.vtk" into this directory
% Run this script


clear all;

% initialize all variables
VTKData.header = cell(5,1);
VTKData.gridType = [];
VTKData.points = [];
VTKData.radius = [];
VTKData.status = [];
VTKData.concentration = [];
VTKData.fibres = [];
VTKData.phenot = [];
VTKData.contpheno = [];

fid = fopen(["first_simulation_cells.1000.vtk"],'r');

% 5 lines of header including DATASET UNSTRUCTURED GRID
VTKData.header{1} = fgets(fid);
if (strcmp(VTKData.header{1}(3:5),'vtk') ~= 1)
    disp('Verify vtk file type and header structure');
    return;
end
VTKData.header{2} = fgets(fid);
VTKData.header{3} = fgets(fid);
VTKData.header{4} = fgets(fid);
VTKData.header{5} = fgets(fid);
VTKData.gridType = VTKData.header{5}(9:end);

% flags for searching a title line to find where the numbers are
flag = 0;

% read next line
str = fgets(fid);
if (strcmp(str(1:6),'POINTS') == 1)
    % find the number of entries
    for i = 7:length(str)
        if (str(i) ~= ' ' && flag == 0)
            flag = 1;
        end
        if (str(i) == ' ' && flag == 1)
            flag = 0;
            numPoints = str2double(str(7:i-1));
            break;
        end
    end
else
    return;
end
% initialize the data structure for points
VTKData.points = zeros(numPoints,3);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.points(i,:) = sscanf(str,'%f',[1 3]);
end
% skip certain lines
for i = 1:4
    skip = fgets(fid);
end
% initialize the data structure for radius
VTKData.radius = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.radius(i,:) = sscanf(str,'%f');
end
% skip certain lines
for i = 1:3
    skip = fgets(fid);
end
% initialize the data structure for status
VTKData.status = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.status(i) = sscanf(str,'%f');
end
% skip certain lines
for i = 1:3
    skip = fgets(fid);
end
% initialize the data structure for concentration
VTKData.concentration = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.concentration(i) = sscanf(str,'%f');
end
% skip certain lines
for i = 1:3
    skip = fgets(fid);
end
% initialize the data structure for fibres
VTKData.fibres = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.fibres(i,:) = sscanf(str,'%f');
end
% skip certain lines
for i = 1:3
    skip = fgets(fid);
end
% initialize the data structure for phenot
VTKData.phenot = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.phenot(i,:) = sscanf(str,'%f');
end
% skip certain lines
for i = 1:3
    skip = fgets(fid);
end
% initialize the data structure for contpheno
VTKData.contpheno = zeros(numPoints,1);
% fill data stucture
for i = 1:numPoints
    str = fgets(fid);
    VTKData.contpheno(i) = sscanf(str,'%f');
end
fclose(fid);

