function varargout = tc_params(varargin)
% TC_PARAMS MATLAB code for tc_params.fig
%      TC_PARAMS, by itself, creates a new TC_PARAMS or raises the existing
%      singleton*.
%
%      H = TC_PARAMS returns the handle to a new TC_PARAMS or the handle to
%      the existing singleton*.
%
%      TC_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TC_PARAMS.M with the given input arguments.
%
%      TC_PARAMS('Property','Value',...) creates a new TC_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tc_params_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tc_params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tc_params

% Last Modified by GUIDE v2.5 27-Jul-2018 16:29:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tc_params_OpeningFcn, ...
                   'gui_OutputFcn',  @tc_params_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
               %List of normalized attenuation coefficient for frequency ranges beginning
%at 0.3MHz - from Fry and Barger 1978
%From 2.1 - 4.0 is extrapolated
c03 = 1;
c04 = 10^-(2/20);
c05 = 10^-(3.5/20);
c06 = 10^-(5/20);
c07 = 10^-(6/20);
c08 = 10^-(7.5/20);
c09 = 10^-(9/20);
c10 = 10^-(10/20);
c11 = 10^-(11.5/20);
c12 = 10^-(13/20);
c13 = 10^-(15/20);
c14 = 10^-(16.5/20);
c15 = 10^-(17.5/20);
c16 = 10^-(18.5/20);
c17 = 10^-(19.5/20);
c18 = 10^-(20.5/20);
c19 = 10^-(22/20);
c20 = 10^-(24/20);
c21 = 10^-(26/20);
c22 = 10^-(28/20);
c23 = 10^-(30/20);
c24 = 10^-(32/20);
c25 = 10^-(34/20);
c26 = 10^-(36/20);
c27 = 10^-(38/20);
c28 = 10^-(40/20);
c29 = 10^-(42/20);
c30 = 10^-(44/20);
c31 = 10^-(46/20);
c32 = 10^-(48/20);
c33 = 10^-(50/20);
c34 = 10^-(52/20);
c35 = 10^-(54/20);
c36 = 10^-(56/20);
c37 = 10^-(58/20);
c38 = 10^-(60/20);
c39 = 10^-(62/20);
c40 = 10^-(66/20);

tc_params.attenuation = [c05 c06 c07 c08 c09 c10 c11 c12 c13 c14 c15 c16...
    c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33...
    c34 c35 c36 c37 c38 c39 c40];
assignin('base','tc_params',tc_params);
               
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% End initialization code - DO NOT EDIT


% --- Executes just before tc_params is made visible.
function tc_params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tc_params (see VARARGIN)

% Choose default command line output for tc_params
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tc_params wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tc_params_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function thick_Callback(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thick as text
%        str2double(get(hObject,'String')) returns contents of thick as a double


% --- Executes during object creation, after setting all properties.
function thick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function speed_Callback(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speed as text
%        str2double(get(hObject,'String')) returns contents of speed as a double


% --- Executes during object creation, after setting all properties.
function speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqs_Callback(hObject, eventdata, handles)
% hObject    handle to freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqs as text
%        str2double(get(hObject,'String')) returns contents of freqs as a double


% --- Executes during object creation, after setting all properties.
function freqs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function calc_atten_Callback(hObject, eventdata, handles)
% hObject    handle to calc_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calc_atten as text
%        str2double(get(hObject,'String')) returns contents of calc_atten as a double


% --- Executes during object creation, after setting all properties.
function calc_atten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calc_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Set_Params.
function Set_Params_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tc_params = evalin('base','tc_params');
tc_params.thickness = str2double(handles.thick.String);
tc_params.speed = str2double(handles.speed.String);
tc_params.freq_range = str2num(handles.freqs.String);


% c05 = 1;
% c06 = 10^-(8/20);
% c07 = 10^-(10/20);
% c08 = 10^-(12/20);
% c09 = 10^-(15/20);
% c10 = 10^-(17/20);
% c11 = 10^-(19/20);
% c12 = 10^-(24/20);
% c13 = 10^-(30/20);
% c14 = 10^-(33/20);
% c15 = 10^-(35/20);
% c16 = 10^-(37/20);
% c17 = 10^-(39/20);
% c18 = 10^-(41/20);
% c19 = 10^-(45/20);
% c20 = 10^-(49/20);
% c21 = 10^-(53/20);
% c22 = 10^-(57/20);
% c23 = 10^-(61/20);
% c24 = 10^-(65/20);
% c25 = 10^-(69/20);
% c26 = 10^-(73/20);
% c27 = 10^-(77/20);
% c28 = 10^-(81/20);
% c29 = 10^-(85/20);
% c30 = 10^-(90/20);
% c31 = 10^-(95/20);
% c32 = 10^-(100/20);
% c33 = 10^-(105/20);
% c34 = 10^-(110/20);
% c35 = 10^-(115/20);
% c36 = 10^-(120/20);
% c37 = 10^-(125/20);
% c38 = 10^-(130/20);
% c39 = 10^-(135/20);
% c40 = 10^-(140/20);
% 
% tc_params.attenuation = [c05 c06 c07 c08 c09 c10 c11 c12 c13 c14 c15 c16...
%     c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33...
%     c34 c35 c36 c37 c38 c39 c40];

atten_idx = round(([tc_params.freq_range(1) tc_params.freq_range(2)])*10-4);
tc_params.freq_divisor = tc_params.attenuation(atten_idx(2))/tc_params.attenuation(atten_idx(1));
% param = evalin('base','param');
% param.tc_params = tc_params;
assignin('base','tc_params',tc_params);
    
    


% --- Executes on button press in calc_attens.
function calc_attens_Callback(hObject, eventdata, handles)
tc_params = evalin('base','tc_params');
tc_params.freq_range = str2num(handles.freqs.String);
atten_idx = round(([tc_params.freq_range(1) tc_params.freq_range(2)])*10-4);
tc_params.freq_divisor = tc_params.attenuation(atten_idx(2))/tc_params.attenuation(atten_idx(1));
handles.calc_atten.String = [num2str(round(tc_params.freq_divisor*100,3,'significant')) '%'];
% hObject    handle to calc_attens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in atten_win.
function atten_win_Callback(hObject, eventdata, handles)
tc_params = evalin('base','tc_params');
tc_params.freq_range = str2num(handles.freqs.String);

atten_idx = round(([tc_params.freq_range(1) tc_params.freq_range(2)])*10-4);
tc_params.freq_divisor = tc_params.attenuation(atten_idx(2))/tc_params.attenuation(atten_idx(1));
a = logspace(tc_params.freq_divisor,1,1000)/10;
x = linspace(tc_params.freq_range(1),tc_params.freq_range(2),1000);
axes(handles.atten_axes)
plot(x,a)
handles.atten_axes.XLabel.String = 'Frequency MHz';
handles.atten_axes.YLabel.String = 'Normalized Attenuation %';
title('Compression Filter Window');
handles.calc_atten.String = [num2str(round(tc_params.freq_divisor*100,3,'significant')) '%'];
% hObject    handle to atten_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
