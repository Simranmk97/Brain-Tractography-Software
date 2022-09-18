%% Sorting the nodes of same tract in indivisual matrices

% The .pts files contain the positions (x,y,z) of each node in 3D space, and each row represents a node.
% The .edge files are organized in a way that the rows in the matrix
% corresponds to the same node (row) in the .pts files.

%% Nodes in the same tract 

% When the nodes are in the same tract, the 2nd entry of one row matches
% the 1st entry in the next row, in the .edge files.

%% Identifying Nodes in different tract - Breaks

% When the nodes are in different tract, the 2nd entry of one row does not match
% the 1st entry in the next row, in the .edge files. These are known as
% "breaks". These breaks imply the end of a preceeding tract and
% beginning of one new tract, in the pathway.

%% Organization Task - Group tracts of the same pathways together.

% (TESTING WITH ONE .EDGE FILE)

clear,clc
% Plan - To make the MRI Data usable, convert them to .m files. 
% loading the m files that were converted from the .pts files.

load('RightHyperDirectPathway.m')

X = RightHyperDirectPathway;

% loading the m files that were convereted from the .edge files.

load('RightHyperDirectPathwayEdgeFile.m')

Y = RightHyperDirectPathwayEdgeFile;

%% Find breaks in the matrices (previously .edge files) using the FIND function 

% Plan - Use find function to find the index no. of the 2nd entry of the
% previous row and the 1st entry of the next row. Store the values in two seperate matrices.
% Use find function the 2nd time to find the index numbers of values when 2nd entry
% does not match with 1st entry of next row.
% Then store the index number in a matrix called breaks.
% The breaks matrix gives us the location of the actual breaks of nodes in
% pathways

index_col = find(Y(:,2)<15059); % matrix 1
index_col2 = find(Y(:,2)<15059) + 1; % matrix 2

breaks = find(Y(:,index_col2) ~= Y(:,index_col)) % not working.












