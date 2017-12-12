function varargout = perchZoneGUI(varargin)
% PERCHZONEGUI MATLAB code for perchZoneGUI.fig
%      PERCHZONEGUI, by itself, creates a new PERCHZONEGUI or raises the existing
%      singleton*.
%
%      H = PERCHZONEGUI returns the handle to a new PERCHZONEGUI or the handle to
%      the existing singleton*.
%
%      PERCHZONEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PERCHZONEGUI.M with the given input arguments.
%
%      PERCHZONEGUI('Property','Value',...) creates a new PERCHZONEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before perchZoneGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to perchZoneGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help perchZoneGUI

% Last Modified by GUIDE v2.5 25-Apr-2017 11:24:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @perchZoneGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @perchZoneGUI_OutputFcn, ...
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


% --- Executes just before perchZoneGUI is made visible.
function perchZoneGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to perchZoneGUI (see VARARGIN)

% Set some global variables
global continueAnalysis
global zoneVertices

continueAnalysis=0;
zoneVertices=[];

% Choose default command line output for perchZoneGUI
handles.output = hObject;

slice=varargin{1};
guititle=varargin{2};
h=imagesc(slice);
colormap gray
set(h.Parent.Parent,'Name',guititle)
handles.h=h;
handles.slice=slice;

handles.vertices=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes perchZoneGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = perchZoneGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.vertices;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global continueAnalysis
global zoneVertices

if ~isfield(handles,'isDrawing')
    handles.isDrawing=0;
    isDrawing=0;
    handles.vertices=[];
    vertices=[];
    handles.whichline=[];
    whichline=[];
else
    isDrawing=handles.isDrawing;
    vertices=handles.vertices;
    whichline=handles.whichline;
end
if isDrawing==0
    % Start new line
    if isempty(whichline)
        lastLine=0;
    else
        lastLine=max(whichline);
    end
    currLine=lastLine+1;
else
    currLine=max(whichline);
end
[currVertex_x,currVertex_y]=ginput(1);
currVertex=[currVertex_x currVertex_y];
vertices=[vertices; currVertex];
whichline=[whichline; currLine];
if isDrawing==1
    prevVertex=vertices(size(vertices,1)-1,:);
    line([prevVertex(1) currVertex(1)],[prevVertex(2) currVertex(2)]);
end

% Update handles structure
isDrawing=1;
handles.isDrawing=isDrawing;
handles.vertices=vertices;
handles.whichline=whichline;
zoneVertices=handles.vertices;
guidata(hObject, handles);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global continueAnalysis

vertices=handles.vertices;
patchHand=patch(vertices(:,1),vertices(:,2),'c');
set(patchHand,'FaceAlpha',0.2);
continueAnalysis=1;

guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global continueAnalysis
global zoneVertices

continueAnalysis=1;
zoneVertices=handles.vertices;

% Hint: delete(hObject) closes the figure
delete(hObject);
