function varargout = beautify(varargin)
% BEAUTIFY MATLAB code for beautify.fig
%      BEAUTIFY, by itself, creates a new BEAUTIFY or raises the existing
%      singleton*.
%
%      H = BEAUTIFY returns the handle to a new BEAUTIFY or the handle to
%      the existing singleton*.
%
%      BEAUTIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEAUTIFY.M with the given input arguments.
%
%      BEAUTIFY('Property','Value',...) creates a new BEAUTIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before beautify_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to beautify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help beautify

% Last Modified by GUIDE v2.5 09-Nov-2018 10:53:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @beautify_OpeningFcn, ...
                   'gui_OutputFcn',  @beautify_OutputFcn, ...
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


% --- Executes just before beautify is made visible.
function beautify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to beautify (see VARARGIN)

% Choose default command line output for beautify
handles.output = hObject;
set(handles.xpt,'String',1);
set(handles.ypt,'String',1);
set(handles.zpt,'String',1);
set(handles.tpt,'String',1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes beautify wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = beautify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function xmin_Callback(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin as text
%        str2double(get(hObject,'String')) returns contents of xmin as a double


% --- Executes during object creation, after setting all properties.
function xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xmax_Callback(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmax as text
%        str2double(get(hObject,'String')) returns contents of xmax as a double


% --- Executes during object creation, after setting all properties.
function xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in max_range.
function max_range_Callback(hObject, eventdata, handles)
% hObject    handle to max_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = evalin('base','X');
S = size(X);
%Sets min ranges to 1
set(handles.xmin,'String',1);
set(handles.ymin,'String',1);
set(handles.zmin,'String',1);
set(handles.tmin,'String',1);
set(handles.cmin,'String',min(X(:)));
%Sets max ranges 
set(handles.xmax,'String',S(1));
set(handles.ymax,'String',S(2));
set(handles.zmax,'String',S(3));
set(handles.tmax,'String',S(4));
set(handles.cmax,'String',max(X(:)));



% Hint: get(hObject,'Value') returns toggle state of max_range


% --- Executes on button press in use_chop.
function use_chop_Callback(hObject, eventdata, handles)
% hObject    handle to use_chop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_chop


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
if handles.use_chop.Value
    X = evalin('base','X_c');
else
    X = evalin('base','Xfilt');
end
assignin('base','X',X);
s = size(X);
% Writes number of points in each dimension
set(handles.xtext,'String',s(1));
set(handles.ytext,'String',s(2));
set(handles.ztext,'String',s(3));
set(handles.ttext,'String',s(4));
% Sets ranges for each dimension
max_range_Callback(hObject, eventdata, handles)


% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plotorig.
function plotorig_Callback(hObject, eventdata, handles)
% hObject    handle to plotorig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
X = evalin('base','X');
end
N = get(handles.plot2d,'Value');
if N == 1 %XZ
    p1 = str2double(get(handles.ypt,'String'));
    p2 = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(x,p1,y,p2))';
elseif N == 2 %YZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(p1,x,y,p2))';
elseif N == 3 %XY
       p1 = str2double(get(handles.zpt,'String'));
    p2 = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    imag = squeeze(X(x,y,p1,p2)');
elseif N == 4 %TX
    p1 = str2double(get(handles.ypt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(x,p1,p2,y));
elseif N == 5 %TY
       p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(p1,x,p2,y));
elseif N == 6 %TZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(p1,p2,x,y));
end
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
      clear h
    h = hotcold;
end
if handles.plotmod.Value
    axes(handles.axes3)
else
    axes(handles.axes1)
end
imagesc(imag)
colormap(h)
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
if isnan(c)
else
    caxis(c)
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function xpt_Callback(hObject, eventdata, handles)
% hObject    handle to xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpt as text
%        str2double(get(hObject,'String')) returns contents of xpt as a double


% --- Executes during object creation, after setting all properties.
function xpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymin_Callback(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymin as text
%        str2double(get(hObject,'String')) returns contents of ymin as a double


% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymax_Callback(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymax as text
%        str2double(get(hObject,'String')) returns contents of ymax as a double


% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypt_Callback(hObject, eventdata, handles)
% hObject    handle to ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypt as text
%        str2double(get(hObject,'String')) returns contents of ypt as a double


% --- Executes during object creation, after setting all properties.
function ypt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmin_Callback(hObject, eventdata, handles)
% hObject    handle to zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmin as text
%        str2double(get(hObject,'String')) returns contents of zmin as a double


% --- Executes during object creation, after setting all properties.
function zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmax_Callback(hObject, eventdata, handles)
% hObject    handle to zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmax as text
%        str2double(get(hObject,'String')) returns contents of zmax as a double


% --- Executes during object creation, after setting all properties.
function zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zpt_Callback(hObject, eventdata, handles)
% hObject    handle to zpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpt as text
%        str2double(get(hObject,'String')) returns contents of zpt as a double


% --- Executes during object creation, after setting all properties.
function zpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmin_Callback(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmin as text
%        str2double(get(hObject,'String')) returns contents of tmin as a double


% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmax_Callback(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmax as text
%        str2double(get(hObject,'String')) returns contents of tmax as a double


% --- Executes during object creation, after setting all properties.
function tmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tpt_Callback(hObject, eventdata, handles)
% hObject    handle to tpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tpt as text
%        str2double(get(hObject,'String')) returns contents of tpt as a double


% --- Executes during object creation, after setting all properties.
function tpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot2d.
function plot2d_Callback(hObject, eventdata, handles)
% hObject    handle to plot2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot2d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot2d


% --- Executes during object creation, after setting all properties.
function plot2d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot3d.
function plot3d_Callback(hObject, eventdata, handles)
% hObject    handle to plot3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot3d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot3d


% --- Executes during object creation, after setting all properties.
function plot3d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmin_Callback(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmin as text
%        str2double(get(hObject,'String')) returns contents of cmin as a double


% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmax_Callback(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmax as text
%        str2double(get(hObject,'String')) returns contents of cmax as a double


% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cmapmenu.
function cmapmenu_Callback(hObject, eventdata, handles)
% hObject    handle to cmapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmapmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmapmenu


% --- Executes during object creation, after setting all properties.
function cmapmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
X = evalin('base','X');
end
N = get(handles.plot3d,'Value');
if N == 1 %XZt
    p = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(x,p,y,:));
elseif N == 2 %YZt
    p = str2double(get(handles.xpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(p,x,y,:));
elseif N == 3 %XYt
   p = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    t = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(x,y,p,:));
elseif N == 4 %XYz
   p = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    t = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(x,y,:,p));
elseif N == 5 %XZy
   p = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    imag = squeeze(X(x,:,y,p));
    imag = permute(imag,[1 3 2]);
elseif N == 6 %YZx
   p = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    imag = squeeze(X(:,x,y,p));
    imag = permute(imag,[2 3 1]);
end
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
      clear h
    h = hotcold;
end
if handles.plotmod.Value
    axes(handles.axes3)
else
    axes(handles.axes1)
end
for i = t
imagesc(imag(:,:,i)')
colormap(h)
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
if isnan(c)
else
    caxis(c)
end
title(i);
drawnow;
 pause(str2double(handles.pausems.String)/1e3);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in holdmods.
function holdmods_Callback(hObject, eventdata, handles)
% hObject    handle to holdmods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of holdmods



function xmean_Callback(hObject, eventdata, handles)
% hObject    handle to xmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmean as text
%        str2double(get(hObject,'String')) returns contents of xmean as a double


% --- Executes during object creation, after setting all properties.
function xmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymean_Callback(hObject, eventdata, handles)
% hObject    handle to ymean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymean as text
%        str2double(get(hObject,'String')) returns contents of ymean as a double


% --- Executes during object creation, after setting all properties.
function ymean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmean_Callback(hObject, eventdata, handles)
% hObject    handle to zmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmean as text
%        str2double(get(hObject,'String')) returns contents of zmean as a double


% --- Executes during object creation, after setting all properties.
function zmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmean_Callback(hObject, eventdata, handles)
% hObject    handle to tmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmean as text
%        str2double(get(hObject,'String')) returns contents of tmean as a double


% --- Executes during object creation, after setting all properties.
function tmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ext_fig.
function ext_fig_Callback(hObject, eventdata, handles)
% hObject    handle to ext_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ext_fig


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = evalin('base','X');
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];

%Xmod = X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2));

box = [str2double(get(handles.xmean,'String')) str2double(get(handles.ymean,'String')) str2double(get(handles.zmean,'String'))];

if handles.allt.Value == 1
    n = t(1):t(2);
else
    n = t(3);
end
for i = n
    Xtemp(:,:,:,i) = imboxfilt3(X(x(1):x(2),y(1):y(2),z(1):z(2),i),box);%imboxfilt3(Xmod(:,:,:,n),box);
end

X(x(1):x(2),y(1):y(2),z(1):z(2),n) = Xtemp(:,:,:,n);

assignin('base','Xnew',X);


% --- Executes on button press in allt.
function allt_Callback(hObject, eventdata, handles)
% hObject    handle to allt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allt


% --- Executes on button press in plotmod.
function plotmod_Callback(hObject, eventdata, handles)
% hObject    handle to plotmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotmod


% --- Executes on button press in clearc.
function clearc_Callback(hObject, eventdata, handles)
% hObject    handle to clearc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cmax,'String',[]);
set(handles.cmin,'String',[]);
% Hint: get(hObject,'Value') returns toggle state of clearc


% --- Executes on button press in loadparams.
function loadparams_Callback(hObject, eventdata, handles)
% hObject    handle to loadparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = evalin('base','selection');
set(handles.xmin,'String',s.x(1));
set(handles.xmax,'String',s.x(2));
set(handles.xpt,'String',s.x(3));
set(handles.ymin,'String',s.y(1));
set(handles.ymax,'String',s.y(2));
set(handles.ypt,'String',s.y(3));
set(handles.zmin,'String',s.z(1));
set(handles.zmax,'String',s.z(2));
set(handles.zpt,'String',s.z(3));
set(handles.tmin,'String',s.t(1));
set(handles.tmax,'String',s.t(2));
set(handles.tpt,'String',s.t(3));
if isnan(s.c)
    s.c = [];
    set(handles.cmin,'String',s.c);
    set(handles.cmax,'String',s.c);
else
set(handles.cmin,'String',s.c(1));
set(handles.cmax,'String',s.c(2));
end



% --- Executes on button press in saveparams.
function saveparams_Callback(hObject, eventdata, handles)
% hObject    handle to saveparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection.x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
selection.y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
selection.z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
selection.t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];
selection.c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
assignin('base','selection',selection);



function pausems_Callback(hObject, eventdata, handles)
% hObject    handle to pausems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pausems as text
%        str2double(get(hObject,'String')) returns contents of pausems as a double


% --- Executes during object creation, after setting all properties.
function pausems_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pausems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = evalin('base','X');
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];
if handles.allt.Value
    Xnew = abs(hilbert(X(x(1):x(2),y(1):y(2),z(1):z(2),:)));
else  
    Xnew = abs(hilbert(X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2))));
    %Xnew = abs(X);
end
assignin('base','Xnew',Xnew);