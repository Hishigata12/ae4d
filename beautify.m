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

% Last Modified by GUIDE v2.5 08-Feb-2019 15:20:43

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
set(handles.xslide,'Min',0);
set(handles.ypt,'String',1);
set(handles.yslide,'Min',0);
set(handles.zpt,'String',1);
set(handles.zslide,'Min',0);
set(handles.tpt,'String',1);
set(handles.tslide,'Min',0);
set(handles.use_chop,'Value',1);
load_Callback(hObject, eventdata, handles);

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
if handles.plotmod.Value == 1
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
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
if length(S) >= 4
    set(handles.tmax,'String',S(4));
else
    set(handles.tmax,'String',1);
end
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
    X = evalin('base','Jrecon');
end
X = real(X);
assignin('base','X',X);
s = size(X);
% Writes number of points in each dimension
set(handles.xtext,'String',s(1));
set(handles.xslide,'Max',s(1))
set(handles.xslide,'SliderStep',[1/s(1) 5/s(1)]);
set(handles.ytext,'String',s(2));
set(handles.yslide,'Max',s(2))
set(handles.yslide,'SliderStep',[1/s(2) 5/s(2)]);
set(handles.ztext,'String',s(3));
set(handles.zslide,'Max',s(3))
set(handles.zslide,'SliderStep',[1/s(3) 5/s(3)]);
if length(s) == 4
    set(handles.ttext,'String',s(4));
    set(handles.tslide,'Max',s(4))
    set(handles.tslide,'SliderStep',[1/s(4) 5/s(4)]);
else
    set(handles.ttext,'String',1);
end


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
if handles.plotmod.Value
% handles.axes3.Children.CData = imag;
else
%     handles.axes1.Children.CData = imag;
end
colormap(h)
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
if isnan(c)
else
    caxis(c);
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


% --- Executes on button press in moviebutton.
function moviebutton_Callback(hObject, eventdata, handles)
% hObject    handle to moviebutton (see GCBO)
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
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];

hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
    clear h
    h = hotcold;
end
if handles.savebox.Value
    figure;
    vidwrite(imag(:,:,t),handles,c,h);
else
    if handles.plotmod.Value
        axes(handles.axes3)
    else
        axes(handles.axes1)
    end
    for i = t
        if i == t(1)
                imagesc(imag(:,:,i)');
        else
            if handles.plotmod.Value
                handles.axes3.Children.CData  = imag(:,:,i)';
            else
                handles.axes1.Children.CData  = imag(:,:,i)';
            end
        end
        colormap(h);
        
        if isnan(c)
        else
            caxis(c);
        end
        title(i);
        %drawnow;
        pause(str2double(handles.pausems.String)/1e3);
    end
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



function med_x_Callback(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_x as text
%        str2double(get(hObject,'String')) returns contents of med_x as a double


% --- Executes during object creation, after setting all properties.
function med_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_y_Callback(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_y as text
%        str2double(get(hObject,'String')) returns contents of med_y as a double


% --- Executes during object creation, after setting all properties.
function med_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_z_Callback(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_z as text
%        str2double(get(hObject,'String')) returns contents of med_z as a double


% --- Executes during object creation, after setting all properties.
function med_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
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
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];

%Xmod = X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2));

box = [str2double(get(handles.med_x,'String')) str2double(get(handles.med_y,'String')) str2double(get(handles.med_z,'String'))];

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

if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];

m1 = [str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n1 = [str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) 0];% get(handles.squarify_box,'Value')];
o1 = [str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
if sum(m1) ~= 3
    m = [1 m1];
else
    m = [0 m1];
end
if sum(n1) ~= 3
    n = [1 n1];
else
    n = [0 n1];
end
if sum(o1) ~= 3
    o = [1 o1];
else
    o = [0 o1];
end
Xtemp = X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2));
X2 = filts3D(Xtemp,m,n,o);
%    X(x(1):round(x(2)*o(1)),y(1):round(y(2)*o(2)),z(1):round(z(2)*o(3)),t(1):t(2)) = X2;
assignin('base','Xnew',X2)


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
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
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



function mean_x_Callback(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_x as text
%        str2double(get(hObject,'String')) returns contents of mean_x as a double


% --- Executes during object creation, after setting all properties.
function mean_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_y_Callback(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_y as text
%        str2double(get(hObject,'String')) returns contents of mean_y as a double


% --- Executes during object creation, after setting all properties.
function mean_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_z_Callback(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_z as text
%        str2double(get(hObject,'String')) returns contents of mean_z as a double


% --- Executes during object creation, after setting all properties.
function mean_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_x_Callback(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_x as text
%        str2double(get(hObject,'String')) returns contents of int_x as a double


% --- Executes during object creation, after setting all properties.
function int_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_y_Callback(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_y as text
%        str2double(get(hObject,'String')) returns contents of int_y as a double


% --- Executes during object creation, after setting all properties.
function int_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_z_Callback(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_z as text
%        str2double(get(hObject,'String')) returns contents of int_z as a double


% --- Executes during object creation, after setting all properties.
function int_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
x = str2num(handles.xmin.String):str2num(handles.xmax.String);
y = str2num(handles.ymin.String):str2num(handles.ymax.String);
z = str2num(handles.zmin.String):str2num(handles.zmax.String);
t = str2num(handles.tmin.String):str2num(handles.tmax.String);
X = X(x,y,z,t);
assignin('base','X_c',X);


% --- Executes on slider movement.
function xslide_Callback(hObject, eventdata, handles)
% hObject    handle to xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.xmin.String) str2double(handles.xmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.xslide.Max = t(2);
handles.xslide.Min = t(1);
set(handles.xslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.xpt,'String',val);
set(handles.xslide,'Value',val);
plotorig_Callback(hObject, eventdata, handles)


%handles.xslide.SliderStep = [1 10];

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function xslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function yslide_Callback(hObject, eventdata, handles)
% hObject    handle to yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.ymin.String) str2double(handles.ymax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.yslide.Max = t(2);
handles.yslide.Min = t(1);
set(handles.yslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.ypt,'String',val);
set(handles.yslide,'Value',val);
plotorig_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function yslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zslide_Callback(hObject, eventdata, handles)
% hObject    handle to zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.zmin.String) str2double(handles.zmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.zslide.Max = t(2);
handles.zslide.Min = t(1);
set(handles.zslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.zpt,'String',val);
set(handles.zslide,'Value',val);
plotorig_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function zslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function tslide_Callback(hObject, eventdata, handles)
% hObject    handle to tslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.tmin.String) str2double(handles.tmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.tslide.Max = t(2);
handles.tslide.Min = t(1);
set(handles.tslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.tpt,'String',val);
set(handles.tslide,'Value',val);
plotorig_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in TodB.
function TodB_Callback(hObject, eventdata, handles)
% hObject    handle to TodB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
S = sign(X);
D = abs(X);
B = 20*log10(D./max(D(:)));
C  = B.*S;
assignin('base','Xnew',C);


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
d = size(X);
if ndims(X) == 3
    if handles.vol.Value
        d = [d 1];
    else
        d = [d(1) 1 d(2) d(3)];
    end
elseif ndims(X) == 2
    d = [d(1) 1 d(2) 1];
end
x = str2num(handles.xmin.String):str2num(handles.xmax.String);
y = str2num(handles.ymin.String):str2num(handles.ymax.String);
z = str2num(handles.zmin.String):str2num(handles.zmax.String);
t = str2num(handles.tmin.String):str2num(handles.tmax.String);
a.x = x;
a.y = y;
a.depth = z;
a.stime = t;
ax_c.x = 1:d(1);
ax_c.y = 1:d(2);
ax_c.depth = 1:d(3);
ax_c.stime = 1:d(4);
assignin('base','ax_c',a);


% --- Executes on button press in vol.
function vol_Callback(hObject, eventdata, handles)
% hObject    handle to vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vol


% --- Executes on button press in savebox.
function savebox_Callback(hObject, eventdata, handles)
% hObject    handle to savebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savebox



function savefolder_Callback(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefolder as text
%        str2double(get(hObject,'String')) returns contents of savefolder as a double


% --- Executes during object creation, after setting all properties.
function savefolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savefigname_Callback(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefigname as text
%        str2double(get(hObject,'String')) returns contents of savefigname as a double


% --- Executes during object creation, after setting all properties.
function savefigname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in thisfold.
function thisfold_Callback(hObject, eventdata, handles)
% hObject    handle to thisfold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.savefolder,'String',pwd);



function framerate_Callback(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate as text
%        str2double(get(hObject,'String')) returns contents of framerate as a double


% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
X = evalin('base','X');
D = size(X);
if length(D) == 2
    D(3:4) = 1;
elseif length(D) == 3
    D(4) = 1;
end
m1 = [str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n1 = [str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) 0];% get(handles.squarify_box,'Value')];
o1 = [str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
if sum(m1) ~= 3
    m = [1 m1];
else
    m = [0 m1];
end
if sum(n1) ~= 3
    n = [1 n1];
else
    n = [0 n1];
end
if sum(o1) ~= 3
    o = [1 o1];
else
    o = [0 o1];
end
n2 = [0 n1];
X2 = filts3D(X,m,n2,o);
S = sign(X2);
X2 = X2.^2;
X2 = X2.*S;
X3 = filts3D(X2,m,n,o);

assignin('base','Xnew',X3);
max_range_Callback(hObject, eventdata, handles)
pushbutton11_Callback(hObject, eventdata, handles)
X = evalin('base','Xnew');
X4 = filts3D(X,m,n2,o);
assignin('base','Xnew',X4);
max_range_Callback(hObject, eventdata, handles)
a = 3;

% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    


% --- Executes on button press in shift.
function shift_Callback(hObject, eventdata, handles)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end

if handles.shiftall.Value
    xs = str2double(handles.xshift.String);
    ys = str2double(handles.yshift.String);
    zs = str2double(handles.zshift.String);
    ts = str2double(handles.tshift.String);
    X = circshift(X,xs,1);
    X = circshift(X,ys,2);
    X = circshift(X,zs,3);
    X = circshift(X,ts,4);
else
    xs = str2double(handles.xshift.String);
    ys = str2double(handles.yshift.String);
    zs = str2double(handles.zshift.String);
    ts = str2double(handles.tshift.String);
    xp = str2double(handles.xpt.String);
    X(xp,:,:,:) = circshift(X(xp,:,:,:),zs,3);
end
assignin('base','Xnew',X);

% hObject    handle to shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tshift_Callback(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tshift as text
%        str2double(get(hObject,'String')) returns contents of tshift as a double


% --- Executes during object creation, after setting all properties.
function tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zshift_Callback(hObject, eventdata, handles)
% hObject    handle to zshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zshift as text
%        str2double(get(hObject,'String')) returns contents of zshift as a double


% --- Executes during object creation, after setting all properties.
function zshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yshift_Callback(hObject, eventdata, handles)
% hObject    handle to yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yshift as text
%        str2double(get(hObject,'String')) returns contents of yshift as a double


% --- Executes during object creation, after setting all properties.
function yshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xshift_Callback(hObject, eventdata, handles)
% hObject    handle to xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xshift as text
%        str2double(get(hObject,'String')) returns contents of xshift as a double


% --- Executes during object creation, after setting all properties.
function xshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shiftall.
function shiftall_Callback(hObject, eventdata, handles)
% hObject    handle to shiftall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shiftall


% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
X = evalin('base','X');
end
xs = str2double(handles.xshift.String);
if xs ~= 0
    xr = -xs:xs;
else
    xr = 0;
end
ys = str2double(handles.yshift.String);
if ys ~= 0
    yr = -ys:ys;
else
    yr = 0;
end
zs = str2double(handles.zshift.String);
ts = str2double(handles.tshift.String);
xp = str2double(handles.xpt.String);
yp = str2double(handles.ypt.String);
tp = str2double(handles.tpt.String);
Y = X;
for i = str2double(handles.tmin.String):str2double(handles.tmax.String)
%     for i = tp
%         figure;
        hold all;
    orig = squeeze(X(xp,yp,:,i));
    for j = xr
        for k = yr
            comp = squeeze(X(xp+j,yp+k,:,i));
%             plot(comp);
            test = corr(orig,comp);
           % disp(test);
            if abs(test) > 0.1  
                cor = xcorr(orig,comp);
%                 plot(cor);
                s = round(length(cor)/2)+find(cor == max(cor));
                if s~=0
%                     disp(s)
                    Y(xp+j,yp+k,:,i) = circshift(comp,s);
                end
            end
        end
    end
    multiWaitbar('Shifting Depths',i/str2double(handles.tmax.String));
end
assignin('base','Xnew',Y);
multiWaitbar('CLOSEALL');


% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
