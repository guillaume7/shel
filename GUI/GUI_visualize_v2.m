%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     SHEL is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with SHEL.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function varargout = visualize_v2(varargin)
% VISUALIZE_V2 M-file for visualize_v2.fig
%      VISUALIZE_V2, by itself, creates a new VISUALIZE_V2 or raises the existing
%      singleton*.
%
%      H = VISUALIZE_V2 returns the handle to a new VISUALIZE_V2 or the handle to
%      the existing singleton*.
%
%      VISUALIZE_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZE_V2.M with the given input arguments.
%
%      VISUALIZE_V2('Property','Value',...) creates a new VISUALIZE_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualize_v2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualize_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visualize_v2

% Last Modified by GUIDE v2.5 13-May-2011 15:39:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualize_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @visualize_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before visualize_v2 is made visible.
function visualize_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualize_v2 (see VARARGIN)

% Choose default command line output for visualize_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualize_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visualize_v2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
    
        s = model_handles;
        s.model_compute_output(handles);
%
        
% --- Executes on button press in printbutton.
function printbutton_Callback(hObject, eventdata, handles)
% hObject    handle to printbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global film;
    global printG;
    global frame;
    global optprint;
    global myfile;
    global saveconf;
    global mydir;
    global myfullfile;
    
    myfile = get(handles.editmyfile,'String');
    myfullfile = myfile;
    mydir = ['data/output/',myfile];
    
    switch optprint        
        case {1, 6, 11} %all
            myfullfile = [myfile,'-all'];        
        case {2, 7, 12} %eta
            myfullfile = [myfile,'-eta'];       
        case {3, 8, 13} %uv
            myfullfile = [myfile,'-uv'];       
        case {4, 9, 14} %left
            myfullfile = [myfile,'-l'];       
        case {5, 10, 15} %right
            myfullfile = [myfile,'-r'];
    end
    
    if exist(mydir,'dir') == 0
        mkdir(mydir);
    end
    
    if film
        printG=true;
        s = model_handles;
        s.model_compute_output(handles);
    else
        g = GUI_common_handles;
        g.GUI_printit(frame, handles, optprint);
        frame = frame + 1;
    end
    printG=false;
    
    %When printing *any* image,
    %*always* save (at the end of simulation)
    %the global variables data for later
    %analysis.
    ExportGlobalQuantities
    
    %Printing an image? Then save the configuration file as well, plz ;)
    %activate the messenger
    saveconf = true; 
    %call the save file callback
    fig2 = GUI_ControlPanel2D;
    handles = guihandles(fig2);
    GUI_ControlPanel2D('savebutton_Callback',fig2,0,handles);
    %defuse the messenger
    saveconf = false;
    GUI_visualize_v2;
    
%%Function ExportGlobalQuantities
function ExportGlobalQuantities
%function ExportGlobalQuantities

    %time-dependent global properties
    global time;
    global volume;
    global iKe;
    global iPe;
    global vtime;
    global iVort;
    global iEnst;
    global iSqStrech;
    global iSqShear;
    global iSqStrain;
    global iOWeiss;
    
    %filename
    global myfile;
    global mydir;

    save([mydir,'/',myfile,'-globaldata']);

% --- Executes during object creation, after setting all properties.
function editmyfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmyfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editmyfile_Callback(hObject, eventdata, handles)
% hObject    handle to editmyfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmyfile as text
%        str2double(get(hObject,'String')) returns contents of editmyfile as a double
    global myfile;  
    value = get(hObject,'String');
    myfile = value;

function setjpeg(val_L, handles)
    set(handles.jpegradio, 'Value', val_L);
    set(handles.epsradio, 'Value', ~val_L);
    set(handles.bmpradio, 'Value', false);

% --- Executes on button press in resetbutton.
function resetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to resetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    s = model_handles;
    s.model_compute_input();
    g = GUI_common_handles;
    g.GUI_plotmodel(handles);

% --- Executes during object creation, after setting all properties.
function EtaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EtaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
global CLimEta;
set(hObject,'Min',-10.);
set(hObject,'Max',2.);
set(hObject,'Value',CLimEta);

% --- Executes on slider movement.
function EtaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to EtaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    global CLimEta;
    value = get(hObject,'Value');
    strin = sprintf('%0.5g',value);
    set(handles.textCLimEta,'String',strin);
    CLimEta = value;
    GUI_plotmodel(handles);
    %modelaxes(handles);

% --- Executes during object creation, after setting all properties.
function velocityslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velocityslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
global CLimV;
set(hObject,'Min',0.);
set(hObject,'Max',100.);
set(hObject,'Value',CLimV);

% --- Executes on slider movement.
function velocityslider_Callback(hObject, eventdata, handles)
% hObject    handle to velocityslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    global CLimV;
    value = get(hObject,'Value');
    strin = sprintf('%0.5g m/s',value);
    set(handles.textCLimV,'String',strin);
    CLimV = value;
    modelaxes(handles);


% --- Executes on button press in filmbox.
function filmbox_Callback(hObject, eventdata, handles)
% hObject    handle to filmbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filmbox
    global film;
    value = get(hObject,'Value');
    film = value;

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
%gets the selected option
    global optprint;
    global film;
    optprint = get(handles.popupmenu2,'Value');
    switch optprint
        case {11, 12, 13, 14, 15}
            set(handles.filmbox,'Value',1)
            film = 1;
    end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in upbutton.
function upbutton_Callback(hObject, eventdata, handles)
% hObject    handle to upbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;

    axes(handles.levelaxes);
    camorbit(0.,2.);
    [az_level el_level] = view;
    el_level = el_level*90;

% --- Executes on button press in downbutton.
function downbutton_Callback(hObject, eventdata, handles)
% hObject    handle to downbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;

    axes(handles.levelaxes);
    camorbit(0.,-2.);
    [az_level el_level] = view;
    el_level = el_level*90;

% --- Executes on button press in leftbutton.
function leftbutton_Callback(hObject, eventdata, handles)
% hObject    handle to leftbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;

    axes(handles.levelaxes);
    camorbit(-2.,0.);
    [az_level el_level] = view;
    el_level = el_level*90;

% --- Executes on button press in rightbutton.
function rightbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rightbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;

    axes(handles.levelaxes);
    camorbit(2.,0.);
    [az_level el_level] = view;
    el_level = el_level*90;
    
% --- Executes on button press in inbutton.
function inbutton_Callback(hObject, eventdata, handles)
% hObject    handle to inbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    camzoom(handles.levelaxes, 1.1);

% --- Executes on button press in outbutton.
function outbutton_Callback(hObject, eventdata, handles)
% hObject    handle to outbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    camzoom(handles.levelaxes, 0.9);


% --- Executes on button press in centerbutton.
function centerbutton_Callback(hObject, eventdata, handles)
% hObject    handle to centerbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global az_level;
    global el_level;
    
    axes(handles.levelaxes);
    [az_level el_level] = view(3);
    a=1;

function writeTextCLims(hobject, eventdata, handles)
%function writeTextCLims(hobject, eventdata, handles)
%
% This function writes the z-scale in the respective textboxes.
%
    global CLimEta;

    writeTextCLim(CLimEta,handles.textCLimEta);
    
function writeTextCLim(value, hobj)    
%function writeTextCLim(value)    
    strin = sprintf('%0.5g',value);
    set(hobj,'String',strin);
    
function showupvalues(hObject, handles)

    writeTextCLims(hObject,0,handles);
    
    %presets the correct view
    view(handles.levelaxes,3);
    

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior T�cnico
%da Universidade T�cnica de Lisboa.
%--------------------------------------------------------------------------
    
