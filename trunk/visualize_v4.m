function varargout = visualize_v4(varargin)
% VISUALIZE_V4 M-file for visualize_v4.fig
%      VISUALIZE_V4, by itself, creates a new VISUALIZE_V4 or raises the existing
%      singleton*.
%
%      H = VISUALIZE_V4 returns the handle to a new VISUALIZE_V4 or the handle to
%      the existing singleton*.
%
%      VISUALIZE_V4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZE_V4.M with the given input arguments.
%
%      VISUALIZE_V4('Property','Value',...) creates a new VISUALIZE_V4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualize_v4_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualize_v4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visualize_v4

% Last Modified by GUIDE v2.5 27-Jan-2010 20:18:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualize_v4_OpeningFcn, ...
                   'gui_OutputFcn',  @visualize_v4_OutputFcn, ...
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


% --- Executes just before visualize_v4 is made visible.
function visualize_v4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualize_v4 (see VARARGIN)

% Choose default command line output for visualize_v4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualize_v4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visualize_v4_OutputFcn(hObject, eventdata, handles)
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
    
        ComputeModel(handles);

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
    mydir = ['images-results/',myfile];
    
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
        ComputeModel(handles);
    else
        printit(frame, handles, optprint);
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
    fig2 = ControlPanel2D;
    handles = guihandles(fig2);
    ControlPanel2D('savebutton_Callback',fig2,0,handles);
    %defuse the messenger
    saveconf = false;
    visualize_v4;
    
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

% --- Executes on button press in button.
function resetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to resetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    initialconditions;
    plotmodel(handles);
    
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
    plotmodel(handles);
    %modelaxes(handles);

function writeTextCLims(hobject, eventdata, handles)
%function writeTextCLims(hobject, eventdata, handles)
%
% This function writes the z-scale in the respective textboxes.
%
    global CLimEta;
    global CLimLeft;
    global CLimRight;

    writeTextCLim(CLimEta,handles.textCLimEta);
    writeTextCLim(CLimLeft,handles.textCLimLeft);
    writeTextCLim(CLimRight,handles.textCLimRight);
    
function writeTextCLim(value, hobj)    
%function writeTextCLim(value)    
    strin = sprintf('%0.5g',value);
    set(hobj,'String',strin);
    
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

% --- Executes on slider movement.
function Leftslider_Callback(hObject, eventdata, handles)
% hObject    handle to Leftslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    global CLimLeft;
    value = get(hObject,'Value');
    strin = sprintf('%0.5g',value);
    set(handles.textCLimLeft,'String',strin);
    CLimLeft = value;
    plotmodel(handles);

% --- Executes during object creation, after setting all properties.
function Leftslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Leftslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
    global CLimLeft;

    usewhitebg = 1;
    if usewhitebg
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end
    
    set(hObject,'Min',-20.);
    set(hObject,'Max',10.);
    set(hObject,'Value',CLimLeft);

% --- Executes on slider movement.
function Rightslider_Callback(hObject, eventdata, handles)
% hObject    handle to Rightslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    global CLimRight;
    value = get(hObject,'Value');
    strin = sprintf('%0.5g',value);
    set(handles.textCLimRight,'String',strin);
    CLimRight = value;
    plotmodel(handles);

% --- Executes during object creation, after setting all properties.
function Rightslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rightslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
    global CLimRight;

    usewhitebg = 1;
    if usewhitebg
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

    set(hObject,'Min',-20.);
    set(hObject,'Max',10.);
    set(hObject,'Value',CLimRight);

% --- Executes on selection change in leftlistbox.
function leftlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to leftlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns leftlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from leftlistbox
    global leftprop;
    oldprop = leftprop;
    leftprop = get(hObject,'Value');
    if resetview(oldprop, leftprop)
        view(handles.Leftaxes, 3)
    end
    plotmodel(handles);

% --- Executes during object creation, after setting all properties.
function leftlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leftlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in rightlistbox.
function rightlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to rightlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns rightlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rightlistbox
    global rightprop;
    oldprop = rightprop;
    rightprop = get(hObject,'Value');
    if resetview(oldprop, rightprop)
        view(handles.Rightaxes, 3)
    end   
    plotmodel(handles);

% --- Executes during object creation, after setting all properties.
function rightlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bool = resetview(oldprop, prop)

    switch oldprop
        %1D or 2D plots (plot, pcolor, contour)
        case {1,3,4,5,6,7,8,22,23,24,25,26,27,28,29,30,31,32,33}
            switch prop
                %3D plots (surf)
                case {2,9,10,11,12,13,14,15,16,17,18,19,20,21}
                    bool = true;
                otherwise
                    bool = false;
            end
        otherwise
            bool = false;
    end
    
%% This function loads the previously saved variables in the visualize4
%% panel.
function showupvalues(hObject, handles)

    global leftprop;
    global rightprop;
    global CLimEta;
    global CLimLeft;
    global CLimRight;
    global film;
    global myfile;
    global optprint;
    
    %No need to plotmodel since they-re plotmodeled right after
    %plotmodel(handles);
    
    set(handles.leftlistbox, 'Value', leftprop);
    set(handles.rightlistbox, 'Value', rightprop);
    set(handles.filmbox, 'Value', film);
    set(handles.EtaSlider, 'Value', CLimEta);
    set(handles.Leftslider, 'Value', CLimLeft);
    set(handles.Rightslider, 'Value', CLimRight);
    if ~optprint
        optprint = 1;
    end
    set(handles.popupmenu2, 'Value', optprint);
    set(handles.editmyfile, 'String', myfile);
    
    writeTextCLims(hObject,0,handles);
    
    %presets the correct view
    view(handles.levelaxes,3);
    setview(handles.Leftaxes, leftprop);
    setview(handles.Rightaxes, rightprop);
    
function setview(Hn, prop)
    switch prop
        case {1,3,4,5,6,7,8,32,33} % 2D
            view(Hn,2);
        case {2,9,10,11,12,13,14,15,16,17,18,19,20,21} % 3D
            view(Hn,3);
    end