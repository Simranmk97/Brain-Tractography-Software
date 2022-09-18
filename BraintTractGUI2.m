function varargout = BraintTractGUI2(varargin)
% BRAINTTRACTGUI2 MATLAB code for BraintTractGUI2.fig
%      BRAINTTRACTGUI2, by itself, creates a new BRAINTTRACTGUI2 or raises the existing
%      singleton*.
%
%      H = BRAINTTRACTGUI2 returns the handle to a new BRAINTTRACTGUI2 or the handle to
%      the existing singleton*.
%
%      BRAINTTRACTGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINTTRACTGUI2.M with the given input arguments.
%
%      BRAINTTRACTGUI2('Property','Value',...) creates a new BRAINTTRACTGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BraintTractGUI2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BraintTractGUI2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BraintTractGUI2

% Last Modified by GUIDE v2.5 25-Apr-2020 13:10:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BraintTractGUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @BraintTractGUI2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end
% --- Executes just before BraintTractGUI2 is made visible.
function BraintTractGUI2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BraintTractGUI2 (see VARARGIN)

% Choose default command line output for BraintTractGUI2
handles.output = hObject;

% data for RHD
    function [break_index_X] = RightHyper(RightHyperDirectPathwayEdgeFile, X1)
disp('Select Right Hyper Direct Pathway edge file ')
RightHyperDirectPathwayEdgeFile = load(uigetfile('*.edge')); 

X = RightHyperDirectPathwayEdgeFile;
for i = 2:length(X(:,1));
    if X(i-1,2) ~= X(i,1);
        breaks2(i) = i;
    else
        breaks2(i) = 0;
    end
end
break_index_X = nonzeros(breaks2);
disp('Select Right Hyper Direct Pathway pts file ')
X1 = load(uigetfile('*.pts')); % RightHyperDirectPathway.pts 
w = 1; 
for i = 1:length(break_index_X);
    w(i+1) = break_index_X(i); 
    RHD.tract{i} = X1(w(i)+i-1:w(i+1)+i-1,:); 
end

end
handles.RHD = RightHyper;

% data for RID
function [break_index_Y] = RightIndirect(RightIndirectPathwayEdgeFile, Y1)
disp('Select Right Indirect Pathway edge file ')
RightIndirectPathwayEdgeFile = load(uigetfile('*.edge')); % uigetfile will allow you to load selected files from the file explorer, just make sure you change the file viewer so that all files are included

Y = RightIndirectPathwayEdgeFile;
for i=2:length(Y(:,1)) % from index 2 to the largest dimension of the 1st col in all rows 
    
    if Y(i-1,2) ~= Y(i,1)
        breaks(i) = Y(i,2); % STORING Y(i,2) IN BREAKS ENSURES THAT BREAKS CONTAINS THE INDEX FOR THE START OF EACH NEW TRACT
    else
        breaks(i) = 0;
    end
    
end
breaks;
break_index_Y = nonzeros(breaks); % This removes zeros from breaks. Also, this break index indicates the beginning of a new node.
disp('Select Right Indirect Pathway pts file ')
Y1 = load(uigetfile('*.pts')); % RightIndirectPathway.pts 

k = 1;
for i = 1:length(break_index_Y)
    k(i+1) = break_index_Y(i);
    RID.tract{i} = Y1(k(i):k(i+1)-1,:); % STORES EACH TRACT FROM Y1 BASED ON START INDICES FOR EACH TRACT STORED IN k
end
end


handles.RID = RightIndirect;

% data for RIC
function [break_index_Z] = RightInternal(RightInternalCapsuleEdgeFile, Z1)
disp('Select Right Internal Capsule edge file ')
RightInternalCapsuleEdgeFile = load(uigetfile('*.edge'));

Z = RightInternalCapsuleEdgeFile;
for i = 2:length(Z(:,1))
    if Z(i-1,2) ~= Z(i,1)
        breaks3(i) = i;
    else
        breaks3(i) = 0;
    end
end

break_index_Z = nonzeros(breaks3);
disp('Select Right Internal Capsule pts file ')
Z1 = load(uigetfile('*.pts')); % RightInternalCapsule.pts 

v = 1;
for i = 1:length(break_index_Z)
    v(i+1) = break_index_Z(i);
    RIC.tract{i} = Z1(v(i)+i-1:v(i+1)+i-1,:);
end
end


handles.RIC = RightInternal;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BraintTractGUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = BraintTractGUI2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end
% --- Executes on button press in RHD_pushButton.
function RHD_pushButton_Callback(hObject, eventdata, handles)
% hObject    handle to RHD_pushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TODO: PLOT RHD DATA
% I HAVE DATA LOADED, WHEN RUNNING IT PULLS NEW FIGURE WINDOW. HOW DO I
% PLOT ON AXES WIDGET IN GUI.

% RHD_1 = handles.RHD
% figure
% for i = 1:length(RHD_1.tract)
%     plot3(RHD_1.tract{i}(:,1),RHD_1.tract{i}(:,2),RHD_1.tract{i}(:,3));
%     hold on
% end
%  
% title('Right Hyper Direct Pathway')
% xlabel('X Position in Voxels')
% ylabel('Y Position in Voxels')
% zlabel('Z Position in Voxels')
end
% --- Executes on button press in RID_pushButton.
function RID_pushButton_Callback(hObject, eventdata, handles)
% hObject    handle to RID_pushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% figure
% for i = 1:length(RID.tract)
%     plot3(RID.tract{i}(:,1),RID.tract{i}(:,2),RID.tract{i}(:,3));
%     hold on
% end
% title('Right Indirect Pathway')
% xlabel('X Position in Voxels')
% ylabel('Y Position in Voxels')
% zlabel('Z Position in Voxels')

end
% --- Executes on button press in RIC_pushButton.
function RIC_pushButton_Callback(hObject, eventdata, handles)
% hObject    handle to RIC_pushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% figure
% for i = 1:length(RIC.tract)
%     plot3(RIC.tract{i}(:,1),RIC.tract{i}(:,2),RIC.tract{i}(:,3));
%     hold on
% end
% 
% title('Right Internal Capsule')
% xlabel('X Position in Voxels')
% ylabel('Y Position in Voxels')
% zlabel('Z Position in Voxels')
end
