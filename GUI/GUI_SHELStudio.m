function varargout = GUI_SHELStudio(varargin)
% GUI_SHELStudio MATLAB code for GUI_SHELStudio.fig
%      GUI_SHELStudio, by itself, creates a new GUI_SHELStudio or raises the existing
%      singleton*.
%
%      H = GUI_SHELStudio returns the handle to a new GUI_SHELStudio or the handle to
%      the existing singleton*.
%
%      GUI_SHELStudio('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SHELStudio.M with the given input arguments.
%
%      GUI_SHELStudio('Property','Value',...) creates a new GUI_SHELStudio or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SHELStudio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SHELStudio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_SHELStudio

% Last Modified by GUIDE v2.5 13-May-2011 12:25:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SHELStudio_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SHELStudio_OutputFcn, ...
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


% --- Executes just before GUI_SHELStudio is made visible.
function GUI_SHELStudio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SHELStudio (see VARARGIN)

% Choose default command line output for GUI_SHELStudio
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_SHELStudio wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SHELStudio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function menu_project_Callback(hObject, eventdata, handles)
% hObject    handle to menu_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_project_item_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_project_item_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global az_level;
    global el_level;
    global npanels;

    %GR: defines the az and el view
    az_level = -37.5;
    el_level = 30;    
    
    %GR : Runs the initial conditions script
    %CreateInitialConditions2D_3;    
    g = GUI_common_handles;
    g.GUI_initializedglobals();
    g.GUI_loadsave('load');

    %GR : showsup the values in the control panel
    %showupvalues(handles);
    g.GUI_setPanels(npanels, handles);
    
% --------------------------------------------------------------------
function menu_project_item_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_project_item_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global saveconf
    global myfile
    global mydir
    global loadbathymetry

    %Was the messenger sent telling us a print is occuring?
    g = GUI_common_handles;
    if saveconf
        g.GUI_loadsave('save',[mydir,'/',myfile]);
    else
        g.GUI_loadsave('save');
    end

% --------------------------------------------------------------------
function menu_project_item_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_project_item_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    GUI_about;

% --------------------------------------------------------------------
function menu_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_momentum_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_momentum_item_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_menu_momentum_item_parameters;

% --------------------------------------------------------------------
function menu_tracers_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_common_item_adsolver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_common_item_adsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracer_item_initialcondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracer_item_initialcondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_common_item_boundarycondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_common_item_boundarycondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracer_item_sourcesink_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracer_item_sourcesink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_common_item_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_common_item_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_momentum_item_adsolver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_adsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_momentum_item_initialcondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_initialcondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_menu_momentum_item_initialconditions;

% --------------------------------------------------------------------
function menu_momentum_item_boundarycondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_boundarycondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_menu_momentum_item_boundaryconditions;

% --------------------------------------------------------------------
function menu_momentum_item_forces_Callback(hObject, eventdata, handles)
% hObject    handle to menu_momentum_item_forces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_menu_momentum_item_forces;

% --------------------------------------------------------------------
function menu_mesh_item_horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mesh_item_horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    GUI_SHELStudio_menu_mesh_item_horizontal;

% --------------------------------------------------------------------
function menu_mesh_item_vertical_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mesh_item_vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_project_item_new_Callback(hObject, eventdata, handles)
% hObject    handle to menu_project_item_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracers_item_common_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracers_item_common (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracers_item_temperature_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracers_item_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracers_item_salinity_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracers_item_salinity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_temperature_item_adsolver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_temperature_item_adsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_temperature_item_initialcondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_temperature_item_initialcondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_temperature_item_boundarycondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_temperature_item_boundarycondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_temperature_item_sourcesink_Callback(hObject, eventdata, handles)
% hObject    handle to menu_temperature_item_sourcesink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_temperature_item_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_temperature_item_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tracers_item_generic_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tracers_item_generic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_salinity_item_adsolver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_salinity_item_adsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_salinity_item_initialcondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_salinity_item_initialcondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_salinity_item_boundarycondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_salinity_item_boundarycondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_salinity_item_sourcesink_Callback(hObject, eventdata, handles)
% hObject    handle to menu_salinity_item_sourcesink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_salinity_item_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_salinity_item_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_generic_item_adsolver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generic_item_adsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_generic_item_initialcondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generic_item_initialcondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_menu_tracers_item_generic_subitem_IC;

% --------------------------------------------------------------------
function menu_generic_item_boundarycondition_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generic_item_boundarycondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_generic_item_sourcesink_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generic_item_sourcesink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_generic_item_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generic_item_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_mesh_item_time_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mesh_item_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    GUI_SHELStudio_menu_mesh_item_time;
    %answer = inputdlg({'time duration (s)','dt (s)'},'Time parameters');
