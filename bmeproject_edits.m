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

RightHyperDirectPathwayEdgeFile = load(uigetfile); % uigetfile will allow you to load selected files from the file explorer

X = RightHyperDirectPathwayEdgeFile;

RightIndirectPathwayEdgeFile = load(uigetfile);

Y = RightIndirectPathwayEdgeFile;

RightInternalCapsuleEdgeFile = load(uigetfile);
Z = RightInternalCapsuleEdgeFile;

%% Finding breaks using for loop

for i=2:length(Y(:,1)) % from index 2 to the largest dimension of the 1st col in all rows 
    
    if Y(i-1,2) ~= Y(i,1)
        breaks(i) = i;
    else
        breaks(i) = 0;
    end
    
end
breaks;
break_index_Y = nonzeros(breaks); % This removes zeros from breaks


%% 

for i = 2:length(X(:,1))
    if X(i-1,2) ~= X(i,1)
        breaks2(i) = i;
    else
        breaks2(i) = 0;
    end
end

break_index_X = nonzeros(breaks2);

%%
for i = 2:length(Z(:,1))
    if Z(i-1,2) ~= Z(i,1)
        breaks3(i) = i;
    else
        breaks3(i) = 0;
    end
end

break_index_Z = nonzeros(breaks3);


%% Organizing Task #2 - Extract all the nodes of the same tract and store into indivisual matrices using the break_index data acquired above.

Y1 = load(uigetfile); % RightHyperDirect

% I've tried different ways using while loop, also with a combo of while and find function. Didn't work. 
% I'm not sure anymore how to store
% the nodes seperately.

i = 1;
k = 1;
blength = length(break_index_Y);
y1 = double(length(break_index_Y),3)
while k < length(break_index_Y)
    if i <= break_index_Y
        y1(i) = Y1(i,:);
        i = i+1;
    else
        k = k+1;
    end
end
% 
% find(Y1(34,1))
% 
% find(Y1(i,1)) <= break_index_Y(k) & find(Y1(i,2))<= break_index_Y(k) & find(Y1(i,3))<= break_index_Y(k);





