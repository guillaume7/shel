%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Foobar is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function varargout = ControlPanel2D(varargin)
% CONTROLPANEL2D M-file for ControlPanel2D.fig
%      CONTROLPANEL2D, by itself, creates a new CONTROLPANEL2D or raises the existing
%      singleton*.
%
%      H = CONTROLPANEL2D returns the handle to a new CONTROLPANEL2D or the handle to
%      the existing singleton*.
%
%      CONTROLPANEL2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTROLPANEL2D.M with the given input arguments.
%
%      CONTROLPANEL2D('Property','Value',...) creates a new CONTROLPANEL2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ControlPanel2D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ControlPanel2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ControlPanel2D

% Last Modified by GUIDE v2.5 24-Aug-2010 11:13:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ControlPanel2D_OpeningFcn, ...
                   'gui_OutputFcn',  @ControlPanel2D_OutputFcn, ...
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

% --- Executes just before ControlPanel2D is made visible.
function ControlPanel2D_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ControlPanel2D (see VARARGIN)

% Choose default command line output for ControlPanel2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ControlPanel2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ControlPanel2D_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function editrho0_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editrho_air_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function editNu_CreateFcn(hObject, ~, ~)
% hObject    handle to editNu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editNu_Callback(hObject, ~, handles)
% hObject    handle to editNu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNu as text
%        str2double(get(hObject,'String')) returns contents of editNu as a double
    global Nu;
    value = str2double(get(hObject, 'String'));
    Nu = value;
    updatestability(handles);
    UpdateCaracteristicNumbers(handles);

% --- Executes during object creation, after setting all properties.
function editg_CreateFcn(hObject, ~, ~)
% hObject    handle to editg (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editg_Callback(hObject, ~, handles)
% hObject    handle to editg (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editg as text
%        str2double(get(hObject,'String')) returns contents of editg as a double
    global g;
    value = str2double(get(hObject, 'String'));
    g = value;
    updatestability(handles);
    UpdateCaracteristicNumbers(handles);

% --- Executes during object creation, after setting all properties.
function editlatitude_CreateFcn(hObject, ~, ~)
% hObject    handle to editlatitude (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editlatitude_Callback(hObject, ~, handles)
% hObject    handle to editlatitude (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlatitude as text
%        str2double(get(hObject,'String')) returns contents of editlatitude as a double
    global f;
    value = str2double(get(hObject, 'String'));
    f = getf(value);
    UpdateCaracteristicNumbers(handles);
    
%% function f_L = getf(lat_L) in degrees
function f_L = getf(lat_L)
%function f_L = getf(lat_L) in degrees

    omega_L = .00007272; %rad/s Mellor Intro to Phys Oceanography
    f_L = 2 * omega_L * sin(lat_L*pi/180.);
    
%% --- Executes during object creation, after setting all properties.
function editlb_CreateFcn(hObject, ~, ~)
% hObject    handle to editlb (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editlb_Callback(hObject, ~, ~)
% hObject    handle to editlb (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlb as text
%        str2double(get(hObject,'String')) returns contents of editlb as a double
    global lb;
    value = str2double(get(hObject, 'String'));
    lb = value;

% --- Executes during object creation, after setting all properties.
function edituwind_CreateFcn(hObject, ~, ~)
% hObject    handle to edituwind (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edituwind_Callback(hObject, ~, ~)
% hObject    handle to edituwind (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edituwind as text
%        str2double(get(hObject,'String')) returns contents of edituwind as a double
    global uwind;
    value = str2double(get(hObject, 'String'));
    uwind = value;

% --- Executes during object creation, after setting all properties.
function editvwind_CreateFcn(hObject, ~, ~)
% hObject    handle to editvwind (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editvwind_Callback(hObject, ~, ~)
% hObject    handle to editvwind (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editvwind as text
%        str2double(get(hObject,'String')) returns contents of editvwind as a double
    global vwind;
    value = str2double(get(hObject, 'String'));
    vwind = value;

% --- Executes during object creation, after setting all properties.
function editkarman_CreateFcn(hObject, ~, ~)
% hObject    handle to editkarman (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editkarman_Callback(hObject, ~, ~)
% hObject    handle to editkarman (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editkarman as text
%        str2double(get(hObject,'String')) returns contents of editkarman as a double
    global karman;
    value = str2double(get(hObject, 'String'));
    karman = value;

% --- Executes during object creation, after setting all properties.
function editduration_CreateFcn(hObject, ~, ~)
% hObject    handle to editduration (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editduration_Callback(hObject, ~, ~)
% hObject    handle to editduration (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editduration as text
%        str2double(get(hObject,'String')) returns contents of editduration as a double
    global duration;
    value = str2double(get(hObject, 'String'));
    duration = value;
    updateTime;

% --- Executes during object creation, after setting all properties.
function editdt_CreateFcn(hObject, ~, ~)
% hObject    handle to editdt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editdt_Callback(hObject, ~, handles)
% hObject    handle to editdt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdt as text
%        str2double(get(hObject,'String')) returns contents of editdt as a double
    global dt;
    
    value = str2double(get(hObject, 'String'));
    dt = value;
    updateTime;
    updatestability(handles);

% --- Executes during object creation, after setting all properties.
function editdx_CreateFcn(hObject, ~, ~)
% hObject    handle to editdx (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editdx_Callback(hObject, ~, handles)
% hObject    handle to editdx (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdx as text
%        str2double(get(hObject,'String')) returns contents of editdx as a double
    global dx;
    value = str2double(get(hObject, 'String'));
    dx = value;
    updatestability(handles);
    UpdateCaracteristicNumbers(handles);

% --- Executes during object creation, after setting all properties.
function editdy_CreateFcn(hObject, ~, ~)
% hObject    handle to editdy (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editdy_Callback(hObject, ~, handles)
% hObject    handle to editdy (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdy as text
%        str2double(get(hObject,'String')) returns contents of editdy as a double
    global dy;
    value = str2double(get(hObject, 'String'));
    dy = value;
    updatestability(handles);

% --- Executes on button press in bathymetrybutton.
function bathymetrybutton_Callback(hObject, ~, handles)
% hObject    handle to bathymetrybutton (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bathymetrybutton
    global loadbathymetry;
    value = get(hObject,'Value');
    loadbathymetry = value;

    update_options_causal_relationships(handles);
    
function de_activate_bathymetry(handles, bool)
%function de_activate_bathymetry(handles, bool)
%
%handles structure with access to all user-data and gui-data
%bool boolean variable. When True, loads bathymetry from file
%and hides manual bathymetry creation. Otherwise, does the opposite.
        if bool
            bath = 'on';
            man = 'off';
        else
            bath = 'off';
            man = 'on';
        end
        set(handles.editbathymetry,'Enable', bath)
        set(handles.edit_d0, 'Enable', man)
        set(handles.editN, 'Enable', man)
        set(handles.editM, 'Enable', man)

% --- Executes during object creation, after setting all properties.
function editbathymetry_CreateFcn(hObject, ~, ~)
% hObject    handle to editbathymetry (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editbathymetry_Callback(hObject, ~, ~)
% hObject    handle to editbathymetry (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbathymetry as text
%        str2double(get(hObject,'String')) returns contents of editbathymetry as a double


% --- Executes during object creation, after setting all properties.
function editN_CreateFcn(hObject, ~, ~)
% hObject    handle to editN (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editN_Callback(hObject, ~, handles)
% hObject    handle to editN (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editN as text
%        str2double(get(hObject,'String')) returns contents of editN as a double
    global N;
    value = str2double(get(hObject, 'String'));
    N = value;
    UpdateCaracteristicNumbers(handles);

% --- Executes during object creation, after setting all properties.
function editM_CreateFcn(hObject, ~, ~)
% hObject    handle to editM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editM_Callback(hObject, ~, handles)
% hObject    handle to editM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editM as text
%        str2double(get(hObject,'String')) returns contents of editM as a double
    global M;
    value = str2double(get(hObject, 'String'));
    M = value;
    UpdateCaracteristicNumbers(handles);

% --- Executes on button press in resetbutton.
function resetbutton_Callback(hObject, ~, handles)
% hObject    handle to resetbutton (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;
    global npanels;

    %GR: defines the az and el view
    az_level = -37.5;
    el_level = 30;    
    
    %GR : Runs the initial conditions script
    %CreateInitialConditions2D_3;
    initializeglobals;
    loadsave('load');

    %GR : showsup the values in the control panel
    showupvalues(handles);
    setPanels(npanels, handles);
            
% --- Executes on button press in coriolisbox.
function coriolisbox_Callback(hObject, ~, handles)
% hObject    handle to coriolisbox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of coriolisbox
    global coriolis;
    value = get(hObject, 'Value');
    coriolis = value;
    UpdateCaracteristicNumbers(handles);
    
    update_options_causal_relationships(handles);

function de_activate_coriolis(handles, bool)
    set(handles.editlatitude, 'Enable', bool)
    %set(handles.geostrophicbox, 'Enable', bool)
    
% --- Executes on button press in pressurebox.
function pressurebox_Callback(hObject, ~, handles)
% hObject    handle to pressurebox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of pressurebox
    global pressure;
    value = get(hObject, 'Value');
    pressure = value;
    
    update_options_causal_relationships(handles);
    
function de_activate_pressure(handles,bool)
    set(handles.editg, 'Enable', bool)
    
% --- Executes on button press in windbox.
function windbox_Callback(hObject, ~, handles)
% hObject    handle to windbox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of windbox
    global wind;
    value = get(hObject, 'Value');
    wind = value;
    
    update_options_causal_relationships(handles);

function de_activate_wind(handles,bool)
    set(handles.editrho_air, 'Enable', bool)
    set(handles.editrho0, 'Enable', bool)
    set(handles.edituwind, 'Enable', bool)
    set(handles.editvwind, 'Enable', bool)

% --- Executes on button press in bottombox.
function bottombox_Callback(hObject, ~, handles)
% hObject    handle to bottombox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of bottombox
    global bottom;
    value = get(hObject, 'Value');
    bottom = value;
    
    update_options_causal_relationships(handles);

function de_activate_bottom(handles,bool)
    set(handles.editlb, 'Enable', bool);
    set(handles.editkarman, 'Enable', bool);

% --- Executes on button press in bumpbox.
function bumpbox_Callback(hObject, ~, handles)
% hObject    handle to bumpbox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bumpbox
    global tc_bump;    
    value = get(hObject,'Value');
    tc_bump = value;
    
    update_options_causal_relationships(handles);
    
% --- is called by bumpbox_callback.
function de_activate_bump(handles, bool)
%function de-activate_bump(handles, bool)
% handles structure with handles and user-data
% bool String with 'on' or 'off' as only possible values .
        set(handles.editBumpX, 'Enable',bool);
        set(handles.editBumpY, 'Enable',bool);
        set(handles.editBumpSx, 'Enable',bool);
        set(handles.editBumpSy, 'Enable',bool);
        set(handles.edit_bump_d0, 'Enable',bool);

% --- Executes on button press in geostrophicbox.
function geostrophicbox_Callback(hObject, ~, handles)
% hObject    handle to geostrophicbox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of geostrophicbox
    global tc_geostrophic;    
    value = get(hObject,'Value');
    tc_geostrophic = value;


% --- Executes during object creation, after setting all properties.
function editoutputdt_CreateFcn(hObject, ~, handles)
% hObject    handle to editoutputdt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editoutputdt_Callback(hObject, ~, handles)
% hObject    handle to editoutputdt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editoutputdt as text
%        str2double(get(hObject,'String')) returns contents of editoutputdt as a double
    global outputdt;   
    outputdt = str2double(get(hObject, 'String'));
    updateTime;

function editBumpX_Callback(hObject, ~, handles)
% hObject    handle to editBumpX (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBumpX as text
%        str2double(get(hObject,'String')) returns contents of editBumpX as a double
    global bump_x0;   
    bump_x0 = str2double(get(hObject, 'String'));


% --- Executes during object creation, after setting all properties.
function editBumpX_CreateFcn(hObject, ~, handles)
% hObject    handle to editBumpX (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editBumpY_Callback(hObject, ~, handles)
% hObject    handle to editBumpY (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBumpY as text
%        str2double(get(hObject,'String')) returns contents of editBumpY as a double
    global bump_y0;   
    bump_y0 = str2double(get(hObject, 'String'));


% --- Executes during object creation, after setting all properties.
function editBumpY_CreateFcn(hObject, ~, handles)
% hObject    handle to editBumpY (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBumpSx_Callback(hObject, ~, handles)
% hObject    handle to editBumpSx (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBumpSx as text
%        str2double(get(hObject,'String')) returns contents of editBumpSx as a double
    global bump_sx;   
    bump_sx = str2double(get(hObject, 'String'));
    updatebumpenergy(handles);


% --- Executes during object creation, after setting all properties.
function editBumpSx_CreateFcn(hObject, ~, handles)
% hObject    handle to editBumpSx (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBumpSy_Callback(hObject, ~, handles)
% hObject    handle to editBumpSy (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBumpSy as text
%        str2double(get(hObject,'String')) returns contents of editBumpSy as a double
    global bump_sy;   
    bump_sy = str2double(get(hObject, 'String'));
    updatebumpenergy(handles);

% --- Executes during object creation, after setting all properties.
function editBumpSy_CreateFcn(hObject, ~, handles)
% hObject    handle to editBumpSy (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxIsla.
function checkboxIsla_Callback(hObject, ~, handles)
% hObject    handle to checkboxIsla (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxIsla
    global tc_isla;
    value = get(hObject, 'Value');
    tc_isla = value;
    
    update_options_causal_relationships(handles);
    
function de_activate_isla(handles, bool)
%function de_activate_isla(handles, bool)
%
%handles structure with user-data and gui-data
%bool 'on' or 'off' only possible values
    set(handles.editIslaX, 'Enable', bool)
    set(handles.editIslaY, 'Enable', bool)
    set(handles.editIslaR, 'Enable', bool)

function editIslaX_Callback(hObject, ~, handles)
% hObject    handle to editIslaX (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIslaX as text
%        str2double(get(hObject,'String')) returns contents of editIslaX as a double
    global isla_x0;   
    isla_x0 = str2double(get(hObject, 'String'));


% --- Executes during object creation, after setting all properties.
function editIslaX_CreateFcn(hObject, ~, handles)
% hObject    handle to editIslaX (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editIslaY_Callback(hObject, ~, handles)
% hObject    handle to editIslaY (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIslaY as text
%        str2double(get(hObject,'String')) returns contents of editIslaY as a double
    global isla_y0;   
    isla_y0 = str2double(get(hObject, 'String'));


% --- Executes during object creation, after setting all properties.
function editIslaY_CreateFcn(hObject, ~, handles)
% hObject    handle to editIslaY (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIslaR_Callback(hObject, ~, handles)
% hObject    handle to editIslaR (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIslaR as text
%        str2double(get(hObject,'String')) returns contents of editIslaR as a double
    global isla_R;   
    isla_R = str2double(get(hObject, 'String'));


% --- Executes during object creation, after setting all properties.
function editIslaR_CreateFcn(hObject, ~, handles)
% hObject    handle to editIslaR (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_d0_Callback(hObject, ~, handles)
% hObject    handle to edit_d0 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d0 as text
%        str2double(get(hObject,'String')) returns contents of edit_d0 as a double
    global d0;   
    d0 = str2double(get(hObject, 'String'));
    updatestability(handles);
    UpdateCaracteristicNumbers(handles);

% --- Executes during object creation, after setting all properties.
function edit_d0_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_d0 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_step.
function checkbox_step_Callback(hObject, ~, handles)
% hObject    handle to checkbox_step (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_step
    global step;
    value = get(hObject, 'Value');
    step = value;
    
    update_options_causal_relationships(handles);
    
function de_activate_step(handles, bool)
        set(handles.edit_d0step, 'Enable', bool);

function edit_d0step_Callback(hObject, ~, handles)
% hObject    handle to edit_d0step (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d0step as text
%        str2double(get(hObject,'String')) returns contents of edit_d0step as a double
    global d0_step;   
    d0_step = str2double(get(hObject, 'String'));
    updatestability(handles);

% --- Executes during object creation, after setting all properties.
function edit_d0step_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_d0step (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton2pan.
function radiobutton2pan_Callback(hObject, ~, handles)
% hObject    handle to radiobutton2pan (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2pan
    global npanels;
    handles = guihandles(gcbo);
    value = get(hObject,'Value');
    if value == true
        npanels = 2;
    end
    setPanels(npanels, handles);

% --- Executes on button press in radiobutton4pan.
function radiobutton4pan_Callback(hObject, ~, handles)
% hObject    handle to radiobutton4pan (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4pan
    global npanels;
    handles = guihandles(gcbo);
    value = get(hObject,'Value');
    if value == true
        npanels = 4;
    end
    setPanels(npanels, handles);

function setPanels(npanels, handles)
%function setPanels(npanels, handles)
%2 --> 2 panels
%4 --> 4 panels
    if npanels == 2
        pan2 = true;
    else
        pan2 = false;
    end
    set(handles.radiobutton2pan, 'Value', pan2);
    set(handles.radiobutton4pan, 'Value', ~pan2);

% --- Executes on button press in pushbuttonShow.
function pushbuttonShow_Callback(hObject, ~, handles)
% hObject    handle to pushbuttonShow (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global npanels;
    global fig1;
    global loadbathymetry;
    global bathymetryfile;

    if ishandle(fig1)
        close(fig1);
    end
    
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid');
    
    %GR: loads the bathymetry file
    if loadbathymetry
        bathymetryfile = get(handles.editbathymetry, 'String');
        if exist(bathymetryfile, 'file') == 0
            [filename, pathname] = uigetfile( ...
                                '*.mat', ...
                                'Pick a matlab bathymetry file' ...
                            );
            bathymetryfile = [pathname,filename];
        end
    end
    
    %GR: openup the visualize panels
    if npanels == 4
        fig1 = visualize_v4;        
        %Write up the z-scale in all the panels
        handles = guihandles(fig1);
        visualize_v4('showupvalues',fig1,handles);        
    else
        fig1 = visualize_v2;
        %Write up the z-scale in all the panels
        handles = guihandles(fig1);
        visualize_v2('showupvalues',fig1,handles);        
    end

    hand = guihandles(fig1);
    
    initialconditions;
    plotmodel(hand);

    
% --- Executes on button press in wallbox.
function wallbox_Callback(hObject, ~, handles)
% hObject    handle to wallbox (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wallbox
    global noslip;
    value = get(hObject, 'Value');
    noslip = value;

% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, ~, handles)
% hObject    handle to savebutton (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global saveconf
    global myfile
    global mydir
    global loadbathymetry

    %Was the messenger sent telling us a print is occuring?
    if saveconf
        loadsave('save',[mydir,'/',myfile]);
    else
        loadsave('save');
    end
    
%% initialize global values
function initializeglobals()
%function initializeglobals()
%
%
global rho0;
global rho_air;
global Nu;
global g;
global f;
global lb;
global karman;
global uwind;
global vwind;
global M;
global N;
global L;
global outputL;     %number of iterations for output
global outputdt;
global dx;
global dy;
global dt;
global duration;
global gama;        %gama coefficient for Asselin-Robert filter
global coriolis;    %use coriolis force?
global bottom;      %use bottom stress?
global wind;        %use wind stress?
global pressure;    %use pressure gradient force?
global noslip;      %use no-slip?

%T-cells
global eta0;

%print to file
global jpeg;
global myfile;
global frame;
global bmp;
global film;
global printG;
global optprint;

%CLim coef A: 
global CLimEta; % m
global CLimLeft; % J
global CLimRight; % J

%Testcases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

%bathymetry
global loadbathymetry;
global bathymetryfile;
global step;
global d0;
global d0_step;

%boundary conditions
global u_closed;
global v_closed;
global u_cyclic;
global v_cyclic;

%bump
global bump_x0;
global bump_y0;
global bump_sx;
global bump_sy;
global bump_d0;

%island
global isla_x0;
global isla_y0;
global isla_R;

%tracer
global tracer;
global K;
global TrXo;
global TrYo;
global TrR;
global Treshold;

rho0 = 1000;
rho_air = 1;
Nu = 0;
K = 0;
g = 9.8;
f = 0.;
lb = 0.;
karman = 0.4;
uwind = 0.;
vwind = 0.;
M = 10;
N = 10;
L = 3;
outputL = 1;     %number of iterations for output
outputdt = 2;
dx = 20;
dy = 20;
dt = 10;
duration = 100;
gama = 0.1;        %gama coefficient for Asselin-Robert filter
coriolis = false;    %use coriolis force?
bottom = false;      %use bottom stress?
wind = false;        %use wind stress?
pressure = true;    %use pressure gradient force?
noslip = false;      %use no-slip?

%T-cells
eta0 = 0;

%print to file
jpeg = true;
myfile = 'null.txt';
frame = true;
bmp = false;
film = false;
printG = false;
optprint = false;

%CLim coef A: [-1eA 1eA]
CLimEta = 0.01; % m
CLimLeft = 1.; % J
CLimRight = 1.; % J

%Testcases
tc_bump = false;
tc_geostrophic = false;
tc_isla = false;
tc_taylor = false;

%bathymetry
loadbathymetry = false;
bathymetryfile = '';
step = false;
d0 = 10.;
d0_step = 1.;

%boundary conditions
u_closed = false;
v_closed = false;
u_cyclic = false;
v_cyclic = false;

%bump
bump_x0 = 0;
bump_y0 = 0;
bump_sx = 0;
bump_sy = 0;
bump_d0 = 0.01;

%island
isla_x0 = 0;
isla_y0 = 0;
isla_R = 0;

%island
tracer = false;
K = 1.;
TrXo = 0.;
TrYo = 0.;
TrR = 0.;
Treshold = 0.001;

%% load and save function
function loadsave(atype, afile)
%function loadsave(type, file)
%
%type - 'load' or 'save'
%Initialize global variables script
%       
%
global rho0;
global rho_air;
global Nu;
global g;
global f;
global lb;
global karman;
global uwind;
global vwind;
global M;
global N;
global L;
global outputL; %number of iterations for output
global outputdt;
global dx;
global dy;
global dt;
global duration;
global gama; %gama coefficient for Asselin-Robert filter
global coriolis; %use coriolis force?
global bottom; %use bottom stress?
global wind; %use wind stress?
global pressure; %use pressure gradient force?
global noslip; %use no-slip?

%T-cells
global eta0;

%print to file
global jpeg;
global myfile;
global frame;
global bmp;
global film;
global printG;
global optprint;

%CLim coef A: [-1eA 1eA]
global CLimEta; % m [-10. 2.]
global CLimV; % m/s [0. 100.]
global CLimLeft; % J [-10. 10.]
global CLimRight; % J [-10. 10.]

%Testcases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

%bathymetry
global loadbathymetry;
global bathymetryfile;
global step;
global d0;
global d0_step;

%boundary conditions
global u_closed;
global v_closed;
global flather;
global radiatelevel;
global radiatetan;
global neumann;

%bump
global bump_x0;
global bump_y0;
global bump_sx;
global bump_sy;
global bump_d0;

%island
global isla_x0;
global isla_y0;
global isla_R;

%tracer
global tracer;
global K;
global TrXo;
global TrYo;
global TrR;
global Treshold;

%Panels outputs & properties
global npanels;
global leftprop;
global rightprop;

%The current bathymetry
global d;

disp(nargchk(1, 2, nargin));

switch atype
    case 'load'
        % Kill initial conditions
        clear globals;
        %Eventually change the uiload with
        %the function`*file = uigetfile* then *load(file)*
        uiload;
    case 'save'
        %Eventually change the uisave with
        %the function`*file = uiputfile* then *save(file)*
        switch nargin 
            case 2
                save(afile);
            otherwise
                uisave;
        end
    otherwise
end

%% Shows up the values
function showupvalues(handles)
%function showupvalues(handles)
% --- Is called by resetbutton ----------------

global rho0;
global rho_air;
global Nu;
global g;
global f;
global lb;
global karman;
global uwind;
global vwind;
global M;
global N;
global L;
global outputL;
global dx;
global dy;
global dt;
global coriolis;
global wind;
global bottom;
global pressure;
global noslip;

%Testcases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

%bump
global bump_x0;
global bump_y0;
global bump_sx;
global bump_sy;
global bump_d0;

%island
global isla_x0;
global isla_y0;
global isla_R;

%bath
global loadbathymetry;
global bathymetryfile;
global step;
global d0;
global d0_step;

%boundary conditions
global u_closed;
global v_closed;
global flather;
global radiatelevel;
global radiatetan;
global neumann;

%tracer
global tracer;
global K;
global TrXo;
global TrYo;
global TrR;
global Treshold;

    updatestability(handles);
    updatebumpenergy(handles);
    
    %GR: Sets the checkboxes
    set(handles.bumpbox, 'Value', tc_bump);
    set(handles.geostrophicbox, 'Value', tc_geostrophic);
    set(handles.checkboxIsla, 'Value', tc_isla);
    set(handles.checkboxTaylorCol, 'Value', tc_taylor);
    
    %GR : Sets the the editable parameters
    set(handles.editrho_air, 'String', num2str(rho_air));
    set(handles.editrho0, 'String', num2str(rho0));
    set(handles.editNu, 'String', num2str(Nu));
    set(handles.editg, 'String', num2str(g));
    set(handles.editlatitude, 'String', num2str(getlatitude(f)));
    set(handles.editlb, 'String', num2str(lb));
    set(handles.editkarman, 'String', num2str(karman));
    set(handles.edituwind, 'String', num2str(uwind));
    set(handles.editvwind, 'String', num2str(vwind));
    set(handles.editM, 'String', num2str(M));
    set(handles.editN, 'String', num2str(N));
    set(handles.editduration, 'String', num2str(L*dt));
    set(handles.editoutputdt, 'String', num2str(outputL*dt));
    set(handles.editdx, 'String', num2str(dx));
    set(handles.editdy, 'String', num2str(dy));
    set(handles.editdt, 'String', num2str(dt));
    set(handles.coriolisbox, 'Value', coriolis);
    set(handles.windbox, 'Value', wind);
    set(handles.bottombox, 'Value', bottom);
    set(handles.pressurebox, 'Value', pressure);
    set(handles.wallbox, 'Value', noslip);
    set(handles.editIslaX, 'String', num2str(isla_x0));
    set(handles.editIslaY, 'String', num2str(isla_y0));
    set(handles.editIslaR, 'String', num2str(isla_R));
    set(handles.editBumpX, 'String', num2str(bump_x0));
    set(handles.editBumpY, 'String', num2str(bump_y0));
    set(handles.edit_bump_d0, 'String', num2str(bump_d0));
    set(handles.editBumpSx, 'String', num2str(bump_sx));
    set(handles.editBumpSy, 'String', num2str(bump_sy));
    set(handles.edit_d0, 'String', num2str(d0));
    set(handles.edit_d0step, 'String', num2str(d0_step));
    set(handles.checkbox_step, 'Value', step);
    set(handles.checkboxUClosed, 'Value', u_closed);
    set(handles.checkboxVClosed, 'Value', v_closed);
    set(handles.bathymetrybutton, 'Value', loadbathymetry);
    set(handles.editbathymetry, 'String', bathymetryfile);
    set(handles.checkboxTracer, 'Value', tracer);
    set(handles.editKtr, 'String', num2str(K));
    set(handles.editTrXo, 'String', num2str(TrXo));
    set(handles.editTrYo, 'String', num2str(TrYo));
    set(handles.editTrR, 'String', num2str(TrR));
    set(handles.editTrThreshold, 'String', num2str(Treshold));    
    set(handles.popupmenuLevelOB,'Value', ...
        radiatelevel * 2 + ~radiatelevel * 1);
    set(handles.popupmenuUOB,'Value', ...
        flather * 3 + neumann * 2 + ~(flather || neumann) * 1);
    set(handles.popupmenuVOB,'Value', ...
        radiatetan * 2 + ~radiatetan * 1)
   
    UpdateCaracteristicNumbers(handles);
    
    update_options_causal_relationships(handles);

%% update stability
function updatestability(handles)
%function updatestability(handles)
%GR : updates the stability status

    global dt;
    global dy;
    global dx;
    global d0;
    global d0_step;
    global g;
    global Nu;
    global tracer;
    global K;

    cfl_U = CFL(dt,dx,dy,max(d0,d0_step),g,Nu);
    
    if tracer
        cfl_T = CFL_T(dt,dx,dy,max(d0,d0_step),g,K);
        set(handles.textTrAC, 'String', num2str( cfl_T(1), '%0.3g'));
        set(handles.textTrB, 'String', num2str( cfl_T(2), '%0.3g'));
        set(handles.textTrUmax, 'String', num2str( cfl_T(3), '%0.3g'));
    end
    
    if cfl_U(1)
        strin = 'Ok';
    else
        strin = 'CFL criterion violation';
    end
     
    %GR : sets the numerical status
    set(handles.statustext, 'String', strin);
        
    %GR : Sets the CFL number
    set(handles.couranttext, 'String', num2str(cfl_U(2),'%0.3g'));
    set(handles.peclettext, 'String', num2str(cfl_U(3),'%0.3g'));

%% function lat = getlatitude(f) in degrees
function lat_L = getlatitude(f_L)
%function lat = getlatitude(f) in degrees

    omega_L = .00007272; %rad/s Mellor Intro to Phys Oceanography
    lat_L = 180. / pi * asin( .5 * f_L / omega_L  );
    
%% cfl criterion
function cfl_T = CFL_T(dt, dx, dy, H, g, K)
    %Pt+dt = A Pi+1 + B Pi + C Pi-1
    % 0 <= A,B,C <= 1
    %CFL computed for an upwind tracer

    %Estimated max flow velocity, based on the wave
    %wake velocity of the flow for ripples at the limit
    %of the hydrostatic approximation.
    h = H / 100; %Wave height at the limit of the hydrostatic appr.
    v_est =  h / 2 * sqrt(g / H);

    Cr = v_est * dt ./ min([dx dy]);
    Dif = K * dt ./ min([dx dy]).^2;
    AC = Cr + Dif;
    B = 1 - 4 * AC; %Instabilities (UP or CD)
    peclet = Cr / 2 / Dif; %Negative oscillations (CD)
    Umax = .25 * dx / dt - K / dx; %Max vel (UP)
    
    cfl_T = [B peclet Umax];

function cfl_L = CFL(dt_L, dx_L, dy_L, H_L, g_L, nu)
%function bool_L = CFL(dt_L, dx_L, dy_L, H_L, g_L, nu)
%
%CFL criterion enumerated in Kantha-Clayson p.282 ~wrong
%for the shallow water equations solver.

 c = sqrt( g_L * max(max(H_L)) );

 transportivity = 2 * c * dt_L * ( dx_L^(-1) + dy_L^(-1) );
 
 %Let us estimate v_max based on the flow speed in the wake of ripples
 %travelling across the domain instead
 v_max =  H_L / 100 * sqrt(g_L / H_L);
 %v_max = max([dx_L dy_L]) / dt_L;
 
 positivity = 0.5 * v_max * max([dx_L dy_L]) ./ nu;
 
 if transportivity < 1 && positivity < 1
     cfl_L = [true, transportivity, positivity];
 else
     cfl_L = [false, transportivity, positivity];
 end

function edit_bump_d0_Callback(hObject, ~, handles)
% hObject    handle to edit_bump_d0 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bump_d0 as text
%        str2double(get(hObject,'String')) returns contents of edit_bump_d0 as a double
    global bump_d0;
    value = str2double(get(hObject, 'string'));
    bump_d0 = value;
    updatebumpenergy(handles);

function updatebumpenergy(handles)
    global bump_d0;
    global bump_sx;
    global bump_sy;
    global rho0;
    global g;
    global d0;
    set(handles.bumpenergytext, 'String', ...
        num2str(.25 * rho0 * g * bump_d0^2 * bump_sx * bump_sy * pi,'%0.3g'));
    set(handles.bumpUotext, 'String', ...
        num2str(.5 * bump_d0 * sqrt(g/d0),'%0.1e'));

% --- Executes during object creation, after setting all properties.
function edit_bump_d0_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_bump_d0 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function UpdateCaracteristicNumbers(handles)
%This function updates the Reynolds and the Rossby numbers

    %read these variables
    global f;
    global g; %does not has an update function
    global Nu;
    global d0;
    global d0_bump;
    global dx;
    global N;
    global M;
    global coriolis;

    %update these variables
    global Ro;
    global Re;
    global Celerity;
    global Uo;
    global RoRadius;
    global CaracTime;
    global HLength;
    global VLength;
    
    %Compute the values
    VLength = d0;
    HLength = max(N * dx, M * dx);
    Celerity = sqrt(g * VLength);
   % Uo = .5 * d0_bump * sqrt( g / (d0_bump + d0));
    Uo = 0.1; %m/s
    %Ro = sqrt(g/VLength)/f;
    Ro = Uo / (f * HLength);
    Re = Celerity * VLength / Nu;
    RoRadius = Celerity / f;
    if coriolis
        CaracTime = max( max( sqrt(VLength/g), 1/f ), HLength/Celerity );
    else
        CaracTime = max( sqrt(VLength/g), HLength/Celerity);
    end
    
    %Write the new values
    set(handles.ReynoldsText, 'String', num2str(Re,'%0.3g'));
    set(handles.CelerityText, 'String', num2str(Celerity,'%0.3g'));
    set(handles.CaracTimeText, 'String', num2str(CaracTime,'%0.3g'));
    set(handles.HSizeText, 'String', num2str(HLength,'%0.3g'));
    set(handles.VSizeText, 'String', num2str(VLength,'%0.3g'));
    if coriolis
        set(handles.RossbyText, 'String', num2str(Ro,'%0.3g'));
        set(handles.RossbyRadiusText, 'String', num2str(RoRadius,'%0.3g'));
    else
        set(handles.RossbyText, 'String', '-');
        set(handles.RossbyRadiusText, 'String', '-');
    end

%% --This function will update all the causal relationships!        
function update_options_causal_relationships(handles)
    
    global u_closed;
    global v_closed;
    global loadbathymetry;
    global coriolis;
    global pressure;
    global wind;
    global bottom;
    global tc_bump;
    global tc_isla;
    global tc_taylor;
    global step;
    global tracer;
    
    %u-cyclic or not? UNCOMMENT THIS IF YOU WANT TO KEEP DEVELOPPING
    %CYCLIC BOUNDARIES
%    if (v_closed)
%        de_activate_v_cyclic(handles, 'off')
%    else
%        de_activate_v_cyclic(handles, 'on')
%    end
    
    %v-cyclic or not? UNCOMMENT THIS IF YOU WANT TO KEEP DEVELOPPING
    %CYCLIC BOUNDARIES
%    if (u_closed)
%        de_activate_u_cyclic(handles, 'off')
%    else
%        de_activate_u_cyclic(handles, 'on')
%    end

    %Do we need OBC?
    if (u_closed && v_closed)
        de_activate_obc(handles, 'off')
    else
        de_activate_obc(handles, 'on')
    end

    %manual bathymetry?
    if (loadbathymetry)
        de_activate_bathymetry(handles, 1);
    else
        de_activate_bathymetry(handles, 0);
    end
    
    %Coriolis?
    if coriolis
        de_activate_coriolis(handles, 'on')
    else
        de_activate_coriolis(handles, 'off')
    end

    %Pressure gradient force?
    if pressure
        de_activate_pressure(handles,'on')
    else
        de_activate_pressure(handles,'off')
    end
    
    %wind?
    if wind
        de_activate_wind(handles,'on');
    else
        de_activate_wind(handles,'off');
    end
    
    %bottom?
    if bottom
        de_activate_bottom(handles,'on');
    else
        de_activate_bottom(handles,'off');
    end
    
    %bump or taylor?
    if (tc_bump || tc_taylor)
        de_activate_bump(handles, 'on')
    else
        de_activate_bump(handles, 'off')
    end
    
    %island?
    if (tc_isla)
        de_activate_isla(handles, 'on');
    else
        de_activate_isla(handles, 'off');
    end
    
    %step?
    if step
        de_activate_step(handles, 'on');
    else
        de_activate_step(handles, 'off');
    end
    
    %tracer?
    if tracer
        de_activate_tracer(handles, 'on');
    else
        de_activate_tracer(handles, 'off');
    end

% --- Executes on button press in checkboxTaylorCol.
function checkboxTaylorCol_Callback(hObject, ~, handles)
% hObject    handle to checkboxTaylorCol (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTaylorCol
    global tc_taylor;
    value = get(hObject, 'Value');
    tc_taylor = value;
    update_options_causal_relationships(handles);
    
%% Tracer code    
function de_activate_tracer(handles, bool)
        set(handles.editKtr, 'Enable', bool);
        set(handles.editTrXo, 'Enable', bool);
        set(handles.editTrYo, 'Enable', bool);
        set(handles.editTrR, 'Enable', bool);
        set(handles.editTrThreshold, 'Enable', bool);

% --- Executes on button press in checkboxTracer.
function checkboxTracer_Callback(hObject, ~, handles)
% hObject    handle to checkboxTracer (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTracer
    global tracer;
    value = get(hObject, 'Value');
    tracer = value;    
    updatestability(handles);
    update_options_causal_relationships(handles);
    
function editKtr_Callback(hObject, ~, handles)
% hObject    handle to editKtr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKtr as text
%        str2double(get(hObject,'String')) returns contents of editKtr as a double
    global K;
    value = str2double(get(hObject, 'String'));
    K = value;    
    updatestability(handles);

% --- Executes during object creation, after setting all properties.
function editKtr_CreateFcn(hObject, ~, ~)
% hObject    handle to editKtr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTrXo_Callback(hObject, ~, ~)
% hObject    handle to editTrXo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrXo as text
%        str2double(get(hObject,'String')) returns contents of editTrXo as a double
    global TrXo;
    value = str2double(get(hObject, 'String'));
    TrXo = value;
    
% --- Executes during object creation, after setting all properties.
function editTrXo_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrXo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTrYo_Callback(hObject, ~, ~)
% hObject    handle to editTrYo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrYo as text
%        str2double(get(hObject,'String')) returns contents of editTrYo as a double
    global TrYo;
    value = str2double(get(hObject, 'String'));
    TrYo = value;    

% --- Executes during object creation, after setting all properties.
function editTrYo_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrYo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTrR_Callback(hObject, ~, ~)
% hObject    handle to editTrR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrR as text
%        str2double(get(hObject,'String')) returns contents of editTrR as a double
    global TrR;
    value = str2double(get(hObject, 'String'));
    TrR = value;    

% --- Executes during object creation, after setting all properties.
function editTrR_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTrThreshold_Callback(hObject, ~, ~)
% hObject    handle to editTrThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrThreshold as text
%        str2double(get(hObject,'String')) returns contents of editTrThreshold as a double
    global Treshold;
    value = str2double(get(hObject, 'String'));
    Treshold = .001 * value; %Theshold is given by user 
                            %in ppt (parts per thousand).

% --- Executes during object creation, after setting all properties.
function editTrThreshold_CreateFcn(hObject, ~, ~)
% hObject    handle to editTrThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuLevelOB.
function popupmenuLevelOB_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuLevelOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuLevelOB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuLevelOB
    global radiatelevel;
    value =  contents{get(hObject,'Value')};
    switch value
        case 'Sommerfeld'
            radiatelevel = true;
        case 'Dirichelet (CI)'
            radiatelevel = false;
    end

% --- Executes during object creation, after setting all properties.
function popupmenuLevelOB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuLevelOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuUOB.
% NORMAL OPEN BOUNDARY CONDITION
function popupmenuUOB_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuUOB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenuUOB
    global flather;
    global neumann;
    value =  get(hObject,'Value');
    switch value
        case 1 %Dirichelet (CI)
            flather = false;
            neumann = false;
            
        case 2 %Neumann (Null)
            flather = false;
            neumann = true;
            
        case 3 %Flather (NVOE)
            flather = true;
            neumann = false;
    end

% --- Executes during object creation, after setting all properties.
function popupmenuUOB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuVOB.
% TANGENT OPEN BOUNDARY CONDITION
function popupmenuVOB_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuVOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuVOB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuVOB
    global radiatetan;
    value =  get(hObject,'Value');
    switch value
        case 1 %Dirichelet (CI)
            radiatetan = false;            
        case 2 %Sommerfeld
            radiatetan = true;            
    end


% --- Executes during object creation, after setting all properties.
function popupmenuVOB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuVOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxUClosed.
function checkboxUClosed_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUClosed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUClosed
    global u_closed;
    value = get(hObject,'Value');
    u_closed = value;
    update_options_causal_relationships(handles);

% --- Executes on button press in checkboxVClosed.
function checkboxVClosed_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxVClosed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxVClosed
    global v_closed;
    value = get(hObject,'Value');
    v_closed = value;
    update_options_causal_relationships(handles);
    
function de_activate_obc(handles, bool)
%function de_activate_bathymetry(handles, bool)
%
%handles structure with access to all user-data and gui-data
%bool boolean variable. When True, loads bathymetry from file
%and hides manual bathymetry creation. Otherwise, does the opposite.
        set(handles.popupmenuLevelOB,'Enable', bool)
        set(handles.popupmenuUOB,'Enable', bool)
        set(handles.popupmenuVOB,'Enable', bool)

% --- Executes on button press in aboutbutton.
function aboutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to aboutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    about;
%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
