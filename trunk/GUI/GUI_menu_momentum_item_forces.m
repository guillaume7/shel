function varargout = GUI_menu_momentum_item_forces(varargin)
% GUI_MENU_MOMENTUM_ITEM_FORCES MATLAB code for GUI_menu_momentum_item_forces.fig
%      GUI_MENU_MOMENTUM_ITEM_FORCES, by itself, creates a new GUI_MENU_MOMENTUM_ITEM_FORCES or raises the existing
%      singleton*.
%
%      H = GUI_MENU_MOMENTUM_ITEM_FORCES returns the handle to a new GUI_MENU_MOMENTUM_ITEM_FORCES or the handle to
%      the existing singleton*.
%
%      GUI_MENU_MOMENTUM_ITEM_FORCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MENU_MOMENTUM_ITEM_FORCES.M with the given input arguments.
%
%      GUI_MENU_MOMENTUM_ITEM_FORCES('Property','Value',...) creates a new GUI_MENU_MOMENTUM_ITEM_FORCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_menu_momentum_item_forces_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_menu_momentum_item_forces_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_menu_momentum_item_forces

% Last Modified by GUIDE v2.5 18-May-2011 16:18:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_menu_momentum_item_forces_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_menu_momentum_item_forces_OutputFcn, ...
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


% --- Executes just before GUI_menu_momentum_item_forces is made visible.
function GUI_menu_momentum_item_forces_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_menu_momentum_item_forces (see VARARGIN)

% Choose default command line output for GUI_menu_momentum_item_forces
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_menu_momentum_item_forces wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_menu_momentum_item_forces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in menu_momentum_item_forces_pushbutton_ok.
function menu_momentum_item_forces_pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_forces_pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    delete(handles.figure1)

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
