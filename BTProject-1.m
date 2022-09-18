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

disp('Select Right Indirect Pathway edge file ')
RightIndirectPathwayEdgeFile = load(uigetfile('*.edge')); % uigetfile will allow you to load selected files from the file explorer, just make sure you change the file viewer so that all files are included

Y = RightIndirectPathwayEdgeFile;

disp('Select Right Hyper Direct Pathway edge file ')
RightHyperDirectPathwayEdgeFile = load(uigetfile('*.edge')); 

X = RightHyperDirectPathwayEdgeFile;

disp('Select Right Internal Capsule edge file ')
RightInternalCapsuleEdgeFile = load(uigetfile('*.edge'));

Z = RightInternalCapsuleEdgeFile;

%% Finding breaks using for loop

% RightIndirectPathway breaks

for i=2:length(Y(:,1)) % from index 2 to the largest dimension of the 1st col in all rows 
    
    if Y(i-1,2) ~= Y(i,1)
        breaks(i) = Y(i,2); % STORING Y(i,2) IN BREAKS ENSURES THAT BREAKS CONTAINS THE INDEX FOR THE START OF EACH NEW TRACT
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
disp('Select Right Indirect Pathway pts file ')
Y1 = load(uigetfile('*.pts')); % RightIndirectPathway.pts 

k = 1;
for i = 1:length(break_index_Y)
    k(i+1) = break_index_Y(i);
    RID.tract{i} = Y1(k(i):k(i+1)-1,:); % STORES EACH TRACT FROM Y1 BASED ON START INDICES FOR EACH TRACT STORED IN k
end

figure
for i = 1:length(RID.tract)
    plot3(RID.tract{i}(:,1),RID.tract{i}(:,2),RID.tract{i}(:,3));
    hold on
end
title('Right Indirect Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')
%% average pathway *use a smoothing function (smooth data funtion)*

% % make the lengths of all RID.tract{1,:} the same
% [R,C] = size(RID.tract(1,:));
% k = 1;
% for i = 1:C %length(RID.tract{k})
%     length = length(RID.tract{c}(:,1));
%     k = k+1;
% end
% max_length = max(length)
%%
Ric = RID.tract; 
avgRIC = zeros([1,3,1]);
figure
s = 1;
    for j = 1:length(Ric{s}(:,1))% row of a given tract
     
        for k = 1:length(Ric) % tract 
           if j<= length(Ric{k}(:,1))
                 RIDX(k) = Ric{k}(j,1);
                 RIDY(k) = Ric{k}(j,2);
                 RIDZ(k) = Ric{k}(j,3);
          
           else 
                 RIDX(k) = 0;
                 RIDY(k) = 0;
                 RIDZ(k) = 0;
           end
           avgRIC(1,1,k) = mean(RIDX(1,:));
           avgRIC(1,2,k) = mean(RIDY(1,:));
           avgRIC(1,3,k) = mean(RIDZ(1,:));
%            plot3(avgRID(k,1),avgRID(k,2),avgRID(k,3),'k.')
%            hold on 
        end

            s = s+1;
    end
B = smooth3(avgRIC);
figure
for i = 1:length(Ric)
    plot3(B(1,1,i),B(1,2,i),B(1,3,i),'k.')
    hold on 
end 
title('Average Right Indirect Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')

%%

% RightHyperDirectPathway.pts
disp('Select Right Hyper Direct Pathway pts file ')
X1 = load(uigetfile('*.pts')); % RightHyperDirectPathway.pts 

w = 1;
for i = 1:length(break_index_X)
    w(i+1) = break_index_X(i); % end 
    RHD.tract{i} = X1(w(i)+i-1:w(i+1)+i-1,:); % omg I solved it
end

figure
for i = 1:length(RHD.tract)
    plot3(RHD.tract{i}(:,1),RHD.tract{i}(:,2),RHD.tract{i}(:,3));
    hold on
end
title('Right Hyper Direct Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')

% TA - I figured out what the issue was. So when there is a break in the .edge file, 
% the second value of the row where the break occurs is the index for the start of 
% the new tract in the .pts file. So, for example, on line 223 in the .edge file a break occurs. 
% The values in row 223 are 225 and 226. In the .pts file, a tract ends on row 225 and a new one begins on 226. 
% I hope that helps, I have attached my edited code with comments on what I changed 
% in all caps to help distinguish from your own comments.
%% average RHD 
Ric = RHD.tract; 
avgRIC = zeros([1,3,1]);
figure
s = 1;
    for j = 1:length(Ric{s}(:,1))% row of a given tract
     
        for k = 1:length(Ric) % tract 
           if j<= length(Ric{k}(:,1))
                 RIDX(k) = Ric{k}(j,1);
                 RIDY(k) = Ric{k}(j,2);
                 RIDZ(k) = Ric{k}(j,3);
          
           else 
                 RIDX(k) = 0;
                 RIDY(k) = 0;
                 RIDZ(k) = 0;
           end
           avgRIC(1,1,k) = mean(RIDX(1,:));
           avgRIC(1,2,k) = mean(RIDY(1,:));
           avgRIC(1,3,k) = mean(RIDZ(1,:));
%            plot3(avgRID(k,1),avgRID(k,2),avgRID(k,3),'k.')
%            hold on 
        end

            s = s+1;
    end
B = smooth3(avgRIC);
figure
for i = 1:length(Ric)
    plot3(B(1,1,i),B(1,2,i),B(1,3,i),'k.')
    hold on 
end 
title('Average Right Hyper Direct Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')
%% 

% RightInternalCapsule.pts
disp('Select Right Internal Capsule pts file ')
Z1 = load(uigetfile('*.pts')); % RightInternalCapsule.pts 

v = 1;
for i = 1:length(break_index_Z)
    v(i+1) = break_index_Z(i);
    RIC.tract{i} = Z1(v(i)+i-1:v(i+1)+i-1,:);
end

figure
for i = 1:length(RIC.tract)
    plot3(RIC.tract{i}(:,1),RIC.tract{i}(:,2),RIC.tract{i}(:,3));
    hold on
end

title('Right Internal Capsule')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')

%% Average RIC
Ric = RIC.tract; 
avgRIC = zeros([1,3,1]);
figure
s = 1;
    for j = 1:length(Ric{s}(:,1))% row of a given tract
     
        for k = 1:length(Ric) % tract 
           if j<= length(Ric{k}(:,1))
                 RIDX(k) = Ric{k}(j,1);
                 RIDY(k) = Ric{k}(j,2);
                 RIDZ(k) = Ric{k}(j,3);
          
           else 
                 RIDX(k) = 0;
                 RIDY(k) = 0;
                 RIDZ(k) = 0;
           end
           avgRIC(1,1,k) = mean(RIDX(1,:));
           avgRIC(1,2,k) = mean(RIDY(1,:));
           avgRIC(1,3,k) = mean(RIDZ(1,:));
%            plot3(avgRID(k,1),avgRID(k,2),avgRID(k,3),'k.')
%            hold on 
        end

            s = s+1;
    end
B = smooth3(avgRIC);
figure
for i = 1:length(Ric)
    plot3(B(1,1,i),B(1,2,i),B(1,3,i),'k.')
    hold on 
end 
title('Average Right Internal Capsule')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')




