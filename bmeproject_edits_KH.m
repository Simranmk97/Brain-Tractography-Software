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
% loading the m files that were convereted from the .edge files.

RightHyperDirectPathwayEdgeFile = load(uigetfile); % uigetfile will allow you to load selected files from the file explorer, just make sure you change the file viewer so that all files are included

X = RightHyperDirectPathwayEdgeFile;

RightIndirectPathwayEdgeFile = load(uigetfile);

Y = RightIndirectPathwayEdgeFile;

RightInternalCapsuleEdgeFile = load(uigetfile);
Z = RightInternalCapsuleEdgeFile;

%% Finding breaks using for loop

% RightIndirectPathway breaks

for i=2:length(Y(:,1)) % from index 2 to the largest dimension of the 1st col in all rows (not including the first row)
    
    if Y(i-1,2) ~= Y(i,1)
        breaks(i) = i;
    else
        breaks(i) = 0;
    end
    
end
breaks;
break_index_Y = nonzeros(breaks); % This removes zeros from breaks. Also, this break index indicates the beginning of a new node.


%% 
% RightHyperDirectPathway breaks

for i = 2:length(X(:,1))
    if X(i-1,2) ~= X(i,1)
        breaks2(i) = i;
    else
        breaks2(i) = 0;
    end
end

break_index_X = nonzeros(breaks2);

%%

%RightInternalCapsulePathway breaks

for i = 2:length(Z(:,1))
    if Z(i-1,2) ~= Z(i,1)
        breaks3(i) = i;
    else
        breaks3(i) = 0;
    end
end

break_index_Z = nonzeros(breaks3);


%% Organizing Task #2 - Extract all the nodes of the same tract and store into indivisual matrices using the break_index data acquired above.

% RightIndirectPathway.pts 

Y1 = load(uigetfile); % RightIndirectPathway.pts 

k = 1;
% the k(i+1) break index is still that of the previous pathway. Need change the k values before the for loop executes. Add
% + 1 to k values.
for i = 1:length(break_index_Y)
    k(i+1) = break_index_Y(i); % (the next row of node (new pathway)) 
    if k(i) == k(i) % this was supposed solve the issue of having the last row of the previous pathway in the new pathway's list. This only worked for the RID.tract(1) and RID.tract(2), but not for the rest. I think break_index_Y nad Y1 indexes are not aligning right.
        k(i) = k(i) + 1;
        RID.tract{i} = Y1(k(i):k(i+1),:); % storing of rows of xyz node from break to break, in a cell structure.
    else
        RID.tract{i} = Y1(k(i):k(i+1),:); 
    end
end
% I still have to figure out how to disconnect the first and last xyz
% points in plots
figure
for i = 1:length(RID.tract)
    plot3(RID.tract{i}(:,1),RID.tract{i}(:,2),RID.tract{i}(:,3)); 
    hold on
end








