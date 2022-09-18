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

for k=2:length(Y(:,1)) % from index 2 to the largest dimension of the 1st col in all rows 
    
    if Y(k-1,2) ~= Y(k,1)
        breaks(k) = Y(k,2); % STORING Y(i,2) IN BREAKS ENSURES THAT BREAKS CONTAINS THE INDEX FOR THE START OF EACH NEW TRACT
    else
        breaks(k) = 0;
    end
    
end
breaks;
break_index_Y = nonzeros(breaks); % This removes zeros from breaks. Also, this break index indicates the beginning of a new node.

%%

% RightHyperDirectPathway breaks

for k = 2:length(X(:,1))
    if X(k-1,2) ~= X(k,1)
        breaks2(k) = k;
    else
        breaks2(k) = 0;
    end
end

break_index_X = nonzeros(breaks2);

%%

%RightInternalCapsulePathway breaks

for k = 2:length(Z(:,1))
    if Z(k-1,2) ~= Z(k,1)
        breaks3(k) = k;
    else
        breaks3(k) = 0;
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
for k = 1:length(RID.tract)
    plot3(RID.tract{k}(:,1),RID.tract{k}(:,2),RID.tract{k}(:,3));
    hold on
end
title('Right Indirect Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')
%% average pathway *use a smoothing function (smooth data funtion)*
Rid = RID.tract; 
avgRID = zeros([1,3,1]);% make the lengths of all RID.tract{1,:} the same
for i = 1:length(Rid(1,:)) % # of tracts
    L(i) = length(Rid{i}); % Length of each tract.
end
max_length = max(L);
% needs to go through each node
    for j = 1:max_length% node of a given tract
     
        for k = 1:length(Rid) % tract # 
           if j<= length(Rid{k})
                 RIDX(k) = Rid{k}(j,1);
                 RIDY(k) = Rid{k}(j,2);
                 RIDZ(k) = Rid{k}(j,3);
          
           else 
                 RIDX(k) = 0;
                 RIDY(k) = 0;
                 RIDZ(k) = 0;
           end
           if RIDX(k)>0
                AVGX(k) = RIDX(k);
           end 
           if RIDY(k)>0
                AVGY(k) = RIDY(k);
           end 
           if RIDZ(k)>0
                AVGZ(k) = RIDZ(k);
           end 
       
           avgRID(j,1,k) = mean(AVGX(1,:));
           avgRID(j,2,k) = mean(AVGY(1,:));
           avgRID(j,3,k) = mean(AVGZ(1,:));

        end

    end
figure
% Using squeeze() makes data for each dimension 1-dimensional so that it is
% compatible with the smoothdata() function
plot3(smoothdata(squeeze(avgRID(1,1,:))),smoothdata(squeeze(avgRID(1,2,:))),smoothdata(squeeze(avgRID(1,3,:))),'k.')
title('Average Right Indirect Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')
%%

% RightHyperDirectPathway.pts
disp('Select Right Hyper Direct Pathway pts file ')
X1 = load(uigetfile('*.pts')); % RightHyperDirectPathway.pts 

w = 1;
for k = 1:length(break_index_X)
    w(k+1) = break_index_X(k); % end 
    RHD.tract{k} = X1(w(k)+k-1:w(k+1)+k-1,:); % omg I solved it
end

figure
for k = 1:length(RHD.tract)
    plot3(RHD.tract{k}(:,1),RHD.tract{k}(:,2),RHD.tract{k}(:,3));
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
Rhd = RHD.tract; 
avgRHD = zeros([1,3,1]);% make the lengths of all RID.tract{1,:} the same
for i = 1:length(Rhd(1,:)) % # of tracts
    L(i) = length(Rhd{i});
end
max_length = max(L);
% needs to go through each node
    for j = 1:max_length% node of a given tract
     
        for k = 1:length(Rhd) % tract 
           if j<= length(Rhd{k}(:,1))
                 RHDX(k) = Rhd{k}(j,1);
                 RHDY(k) = Rhd{k}(j,2);
                 RHDZ(k) = Rhd{k}(j,3);
          
           else 
                 RHDX(k) = 0;
                 RHDY(k) = 0;
                 RHDZ(k) = 0;
           end
           if RHDX(k)>0
                AVGX(k) = RHDX(k);
           end 
           if RHDY(k)>0
                AVGY(k) = RHDY(k);
           end 
           if RHDZ(k)>0
                AVGZ(k) = RHDZ(k);
           end 
       
           avgRHD(1,1,k) = mean(AVGX(1,:));
           avgRHD(1,2,k) = mean(AVGY(1,:));
           avgRHD(1,3,k) = mean(AVGZ(1,:));
        end

    end
figure
plot3(smoothdata(squeeze(avgRHD(1,1,:))),smoothdata(squeeze(avgRHD(1,2,:))),smoothdata(squeeze(avgRHD(1,3,:))),'k.')
title('Average Right Hyper Direct Pathway')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')
%% 

% RightInternalCapsule.pts
disp('Select Right Internal Capsule pts file ')
Z1 = load(uigetfile('*.pts')); % RightInternalCapsule.pts 

v = 1;
for k = 1:length(break_index_Z)
    v(k+1) = break_index_Z(k);
    RIC.tract{k} = Z1(v(k)+k-1:v(k+1)+k-1,:);
end

figure
for k = 1:length(RIC.tract)
    plot3(RIC.tract{k}(:,1),RIC.tract{k}(:,2),RIC.tract{k}(:,3));
    hold on
end

title('Right Internal Capsule')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')

%% Average RIC
Ric = RIC.tract; 
avgRIC = zeros([1,3,1]);% make the lengths of all RID.tract{1,:} the same
for i = 1:length(Ric(1,:)) % # of tracts
    L(i) = length(Ric{i});
end
max_length = max(L);
% needs to go through each node
    for j = 1:max_length% node of a given tract
     
        for k = 1:length(Ric) % tract 
           if j<= length(Ric{k}(:,1))
                 RICX(k) = Ric{k}(j,1);
                 RICY(k) = Ric{k}(j,2);
                 RICZ(k) = Ric{k}(j,3);
          
           else 
                 RICX(k) = 0;
                 RICY(k) = 0;
                 RICZ(k) = 0;
           end
           if RICX(k)>0
                AVGX(k) = RICX(k);
           end 
           if RICY(k)>0
                AVGY(k) = RICY(k);
           end 
           if RICZ(k)>0
                AVGZ(k) = RICZ(k);
           end 
       
           avgRIC(1,1,k) = mean(AVGX(1,:));
           avgRIC(1,2,k) = mean(AVGY(1,:));
           avgRIC(1,3,k) = mean(AVGZ(1,:));
        end 
    end
figure
plot3(smoothdata(squeeze(avgRIC(1,1,:))),smoothdata(squeeze(avgRIC(1,2,:))),smoothdata(squeeze(avgRIC(1,3,:))),'k.')
title('Average Right Internal Capsule')
xlabel('X Position in Voxels')
ylabel('Y Position in Voxels')
zlabel('Z Position in Voxels')




