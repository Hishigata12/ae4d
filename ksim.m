function varargout = ksim(varargin)
% KSIM MATLAB code for ksim.fig
%      KSIM, by itself, creates a new KSIM or raises the existing
%      singleton*.
%
%      H = KSIM returns the handle to a new KSIM or the handle to
%      the existing singleton*.
%
%      KSIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KSIM.M with the given input arguments.
%
%      KSIM('Property','Value',...) creates a new KSIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ksim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ksim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ksim

% Last Modified by GUIDE v2.5 14-Jun-2019 20:12:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ksim_OpeningFcn, ...
                   'gui_OutputFcn',  @ksim_OutputFcn, ...
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


% --- Executes just before ksim is made visible.
function ksim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ksim (see VARARGIN)

% Choose default command line output for ksim
handles.output = hObject;
if ismember('sensor_data',evalin('base','who'))
    sensor_data = evalin('base','sensor_data');
    set(handles.data_menu,'String',fieldnames(sensor_data));
    if isstruct(sensor_data)
        sensor_data = sensor_data.p;
    end
    update_sensormenu(hObject, eventdata, handles, sensor_data);
    if ismember('kgrid',evalin('base','who'))
        kgrid = evalin('base','kgrid');
        update_text(hObject, eventdata, handles, sensor_data, kgrid);
    end
end
set(handles.sensor_xpts,'Visible','off');
set(handles.sensor_xpts_text,'Visible','off');
set(handles.sensor_ypts,'Visible','off');
set(handles.sensor_ypts_text,'Visible','off');
set(handles.source_menu1,'Value',2);
set(handles.source_menu2,'Value',2);
set(handles.source_menu3,'Value',2);
set(handles.medium_layers,'Visible','off')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ksim wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function update_text(hObject, eventdata, handles, sensor_data, kgrid);
sensor_size = size(sensor_data);
new_size = [num2str(sensor_size(1)), ' Sensors ', num2str(sensor_size(2)),' Points'];
set(handles.time_sensor_text,'String',new_size);

dt = kgrid.dt*1e6; %dt in [us]
set(handles.dt_text,'String',[num2str(dt), ' usec/pt']);

function update_sensormenu(hObject, eventdata, handles, sensor_data);
S = linspace(1,size(sensor_data,1),size(sensor_data,1)); %Sensor numbers
set(handles.sensormenu,'String',S);


% --- Outputs from this function are returned to the command line.
function varargout = ksim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in sensormenu.
function sensormenu_Callback(hObject, eventdata, handles)
% hObject    handle to sensormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sensormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sensormenu


% --- Executes during object creation, after setting all properties.
function sensormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_sensor.
function show_sensor_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sensor_data = evalin('base','sensor_data');
val = handles.data_menu.String{handles.data_menu.Value};

%Gets vars depending on parameter
switch val
    case 'p'
        sensor_data = sensor_data.p;
     %   clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];%[-1 1];
        ylab = 'Pressure';
        ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Pressure'
        plot1 = 0;
    case 'ux'
        sensor_data = sensor_data.ux;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
         ylab = 'Velocity X';
          ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'X Velocity';
         plot1 = 0;
    case 'uy'
        sensor_data = sensor_data.uy;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
         ylab = 'Velocity Y';
           ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Y Velocity'
         plot1 = 0;
    case 'p_final'
        sensor_data = sensor_data.p_final;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Pressure';
         plot1 = 1;
    case 'ux_final'
        sensor_data = sensor_data.ux_final;
       % clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Velocity X';
         plot1 = 1;
    case 'uy_final'
        sensor_data = sensor_data.uy_final;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Velocity Y';
         plot1 = 1;
end

% Gets CLims
if isempty(get(handles.display_clims,'String'))
    clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
else
    clims = str2num(handles.display_clims.String);
end

%Plots on Axes 1
axes(handles.axes1)
if plot1 == 0
imagesc(kgrid.t_array,1:size(sensor_data,1),sensor_data, clims);
elseif plot1 == 1
    imagesc(sensor_data,clims);
end
colormap('hotcold');
ylabel(ylab1);    
xlabel(xlab1);
title(tits1);
colorbar;

%Plots on axes 2
if plot1 == 0
axes(handles.axes2)
Sn = get(handles.sensormenu,'Value');
plot(kgrid.t_array,sensor_data(Sn,:))
ylabel(ylab);
xlabel('Time');
end

%Plots on axes 3
if plot1 == 0
sensor = evalin('base','sensor');
axes(handles.axes3);
hold off
scatter(sensor.mask(1,:),sensor.mask(2,:),'k');
hold all;
shift = find(sensor.mask(1,:) == min(sensor.mask(1,:)));
sensor.mask = circshift(sensor.mask,-shift+1,2);
scatter(sensor.mask(1, Sn),sensor.mask(2,Sn),'r','filled');
% scatter(sensor.mask(1, Sn + shift - 1),sensor.mask(2,Sn + shift - 1),'r','filled');
xsize = [min(kgrid.x_vec),max(kgrid.x_vec)];
ysize = [min(kgrid.y_vec),max(kgrid.y_vec)];
xlim(xsize); ylim(ysize);
colorbar;
end


% hObject    handle to show_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sensor.
function sensor_Callback(hObject, eventdata, handles)
sensor = evalin('base','sensor');
kgrid = evalin('base','kgrid');
axes(handles.axes2);
scatter(sensor.mask(1,:),sensor.mask(2,:));
xsize = [min(kgrid.x_vec),max(kgrid.x_vec)];
ysize = [min(kgrid.y_vec),max(kgrid.y_vec)];
xlim(xsize); ylim(ysize);

% hObject    handle to sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in source.
function source_Callback(hObject, eventdata, handles)
source = evalin('base','source');
kgrid = evalin('base','kgrid');
if ismember('p_mask',fieldnames(source))
    source.p0 = source.p0+source.p_mask;
end
xsize = [min(kgrid.x_vec),max(kgrid.x_vec)];
ysize = [min(kgrid.y_vec),max(kgrid.y_vec)];
axes(handles.axes2);
imagesc(kgrid.x_vec,kgrid.y_vec,source.p0);
xlim(xsize); ylim(ysize);
ylabel('X Position')
xlabel('Y Position')
title('Source')
if ismember('p',fieldnames(source))
axes(handles.axes3)
hold off;
smenu = handles.source_menu1.String{handles.source_menu1.Value};
trans = handles.trans_menu.String{handles.trans_menu.Value};
if strcmp(smenu,'Transducer');
    switch trans
        case 'P4-2','P4-1'
            axes(handles.axes3)
             shift = linspace(0,1,size(source.p,1));
           img = figure;
            for i = 1:size(source.p,1)
               pos = [0 1-shift(i) 1 1/size(source.p,1)];
            subplot('Position',pos)
           img = plot(kgrid.t_array(1:size(source.p,2)),source.p(i,:));
           xticks([]);
           yticks([]);
           set(gca,'Color','none')
           xlabel(''); ylabel('');
            end
        case 'H235','H247'
        case '1MHz'
            plot(kgrid.t_array*1e6,source.p)
            xlabel('Time [us]')
            ylabel('Magnitude')
            title('Excitation')
    end
else
plot(kgrid.t_array*1e6,source.p)
xlabel('Time [us]')
ylabel('Magnitude')
title('Excitation')
end

end
% hObject    handle to source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
%import and convert variables from GUI
medium = evalin('base','medium');
source = evalin('base','source');
sensor = evalin('base','sensor');
kgrid = evalin('base','kgrid');
rec = boolean(get(handles.record_movie,'Value'));
mesh = boolean(get(handles.sim_mesh,'Value'));
mask1 = get(handles.sim_mask,'Value');
% moviename = get(handles.rec_name,'String');
moviename = 'poopooman'
if mask1
    mask = 'on';
else
    mask = 'off';
end
pml = boolean(get(handles.sim_pml,'Value'));
scale = str2num(get(handles.sim_scale,'String'));
freq = str2double(get(handles.sim_freq,'String'));
fps = str2double(get(handles.sim_fps,'String'));

%Prep vars for use
fdir = get(handles.record_dir,'String');
fname = get(handles.record_filename,'String');
fname2 = fullfile(fdir,fname);
rec_args = {'RecordMovie',rec,'MovieName',fname2,'PlotScale',scale,'PlotFreq',freq,...
    'MovieProfile','MPEG-4','MovieArgs',{'FrameRate',fps}};
mesh_args = {'MeshPlot', mesh};
mask_args = {'DisplayMask',mask};
pml_args = {'PlotPML',pml}
% input_args = {'PlotScale',scale,'PlotFreq',freq};
input_args = [rec_args mesh_args, mask_args, pml_args];
kgrid.dt = str2double(get(handles.dt_sim,'String'))*1e-6;
if handles.sim_3D.Value
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
else
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end
axes(handles.axes1);
assignin('base','sensor_data',sensor_data);
set(handles.data_menu,'String',fieldnames(sensor_data));

if isstruct(sensor_data)
    sensor_data = sensor_data.p;
end

update_sensormenu(hObject, eventdata, handles, sensor_data);

update_text(hObject, eventdata, handles, sensor_data, kgrid);

% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tp_Callback(hObject, eventdata, handles)
% hObject    handle to tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tp as text
%        str2double(get(hObject,'String')) returns contents of tp as a double


% --- Executes during object creation, after setting all properties.
function tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_tp.
function show_tp_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sensor_data = evalin('base','sensor_data');
if isstruct(sensor_data)
    sensor_data = sensor_data.p;
end
axes(handles.axes1)
imagesc(kgrid.t_array.*1e6,1:size(sensor_data,1),sensor_data, [-1, 1]);
colormap('hotcold');
ylabel('Sensor Number');
xlabel('Time Step');
colorbar;


Tn = str2num(get(handles.tp,'String'));

for i = 1:length(Tn)
    Tp = Tn(i);
    axes(handles.axes2)
    if i == 1
plot(1:size(sensor_data,1),sensor_data(:,Tp))
    else
        handles.axes2.Children.YData = sensor_data(:,Tp);
    end
    t = num2str(kgrid.t_array(Tp).*1e6);
    title(['tp = ', t, ' us']);
ylabel('Pressure');
xlabel('Sensor');

sensor = evalin('base','sensor');
axes(handles.axes3);
hold off
color_vals = (sensor_data(:,Tp)-min(sensor_data(:)))./max(sensor_data(:));
color_vals = sensor_data(:,Tp);
c = [color_vals, zeros(length(color_vals),2)];
if i ==1 
scatter(sensor.mask(1,:),sensor.mask(2,:),[],color_vals);
colormap('hotcold');
else 
    handles.axes3.Children.CData = color_vals;
end
hold all;
shift = find(sensor.mask(1,:) == min(sensor.mask(1,:)));
% scatter(sensor.mask(1,shift - 1),sensor.mask(2, shift - 1),'r','filled');
xsize = [min(kgrid.x_vec),max(kgrid.x_vec)];
ysize = [min(kgrid.y_vec),max(kgrid.y_vec)];
xlim(xsize); ylim(ysize);

end

% hObject    handle to show_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function dt_sim_Callback(hObject, eventdata, handles)
% hObject    handle to dt_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_sim as text
%        str2double(get(hObject,'String')) returns contents of dt_sim as a double


% --- Executes during object creation, after setting all properties.
function dt_sim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dz_Callback(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dz as text
%        str2double(get(hObject,'String')) returns contents of dz as a double


% --- Executes during object creation, after setting all properties.
function dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nx_Callback(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nx as text
%        str2double(get(hObject,'String')) returns contents of Nx as a double


% --- Executes during object creation, after setting all properties.
function Nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ny_Callback(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny as text
%        str2double(get(hObject,'String')) returns contents of Ny as a double


% --- Executes during object creation, after setting all properties.
function Ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nz_Callback(hObject, eventdata, handles)
% hObject    handle to Nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nz as text
%        str2double(get(hObject,'String')) returns contents of Nz as a double


% --- Executes during object creation, after setting all properties.
function Nz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nt_Callback(hObject, eventdata, handles)
% hObject    handle to Nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nt as text
%        str2double(get(hObject,'String')) returns contents of Nt as a double


% --- Executes during object creation, after setting all properties.
function Nt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in make_grid.
function make_grid_Callback(hObject, eventdata, handles)
Nx = str2num(get(handles.Nx,'String'));
Ny = str2num(get(handles.Ny,'String'));
Nz = str2num(get(handles.Nz,'String'));
Nt = str2num(get(handles.Nt,'String'));
dx = str2num(get(handles.dx,'String'))*1e-3;
dy = str2num(get(handles.dy,'String'))*1e-3;
dz = str2num(get(handles.dz,'String'))*1e-3;
dt = str2num(get(handles.dt_sim,'String'))*1e-6;
if Nz == 0
    if Ny == 0
        kWaveGrid(Nx,dx);
    else
        kgrid = kWaveGrid(Nx,dx,Ny,dy);
    end
else
    kgrid = kWaveGrid(Nx,dx,Ny,dy,Nz,dz);
end
kgrid.Nt = Nt;
    kgrid.dt = dt;
kgrid.t_array = 0:dt:Nt*dt;
assignin('base','kgrid',kgrid);
% hObject    handle to make_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_medium.
%Input 0 for full range of xN, yN

function make_medium_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
mode = handles.medium_menu.String{handles.medium_menu.Value};
switch mode
    case 'Homogenous'
        sound_speed = str2double(handles.medium_ss.String)*1e3;
        alpha_coeff = str2double(handles.medium_alpha_coeff.String);
        alpha_power = str2double(handles.medium_alpha_power.String);
        density = str2double(handles.medium_density.String)*1e3;
        medium.sound_speed = sound_speed; %m/s
        medium.alpha_coeff = alpha_coeff; %dB/(MHz^y cm)
        medium.alpha_power = alpha_power;
        medium.density = density;
    case 'Heterogenous'
        BG = ones(kgrid.Nx,kgrid.Ny);
        sound_speed = str2double(handles.medium_ss.String)*1e3*BG;
        alpha_coeff = str2double(handles.medium_alpha_coeff.String)*BG;
        alpha_power = str2double(handles.medium_alpha_power.String);
        density = str2double(handles.medium_density.String)*1e3*BG;
        medium.sound_speed = sound_speed; %m/s
        medium.alpha_coeff = alpha_coeff; %dB/(MHz^y cm)
        medium.alpha_power = alpha_power;
        medium.density = density;
        %Layer 1
        if handles.ly1_on.Value
            xr1 = str2double(handles.ly1_xN.String);
            yr1 = str2double(handles.ly1_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly1_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly1_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly1_ac.String)*rng1;
            x01 = str2double(handles.ly1_x0.String);
            y01 = str2double(handles.ly1_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
        %Layer 2
        if handles.ly2_on.Value
            xr1 = str2double(handles.ly2_xN.String);
            yr1 = str2double(handles.ly2_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly2_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly2_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly2_ac.String)*rng1;
            x01 = str2double(handles.ly2_x0.String);
            y01 = str2double(handles.ly2_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
        %Layer 3
        if handles.ly3_on.Value
            xr1 = str2double(handles.ly3_xN.String);
            yr1 = str2double(handles.ly3_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly3_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly3_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly3_ac.String)*rng1;
            x01 = str2double(handles.ly3_x0.String);
            y01 = str2double(handles.ly3_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
end

% f_max = medium.sound_speed/2/kgrid.dx; %maximum frequency allowed

assignin('base','medium',medium);
% hObject    handle to make_medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_source.
function make_source_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sourcemenu1 = handles.source_menu1.String{handles.source_menu1.Value};
switch sourcemenu1
    case 'none'
        source_1 = zeros(kgrid.Nx,kgrid.Ny);
    case 'Disk'
        disc_magnitude = str2num(get(handles.source_mag1,'String'));
        disc_x_pos = str2num(get(handles.source_xpos1,'String'));
        disc_y_pos = str2num(get(handles.source_ypos1,'String'));
        disc_rad = str2num(get(handles.source_rad1,'String'));
        source_1 = disc_magnitude * makeDisc(kgrid.Nx, kgrid.Ny, disc_x_pos, disc_y_pos, disc_rad);
    case 'Time-Varying Point'
        %delay, duration and period in micro seconds
        %xpos and ypos in grid points
        magnitude = str2num(get(handles.source_mag1,'String'));
        x_pos = str2num(get(handles.source_xpos1,'String'));
        y_pos = str2num(get(handles.source_ypos1,'String'));
        if numel(x_pos) ~= numel(y_pos)
            errordlg('Xpos and Ypos must be same length')
        end
        delay = str2num(get(handles.source_delay1,'String'));
        duration = str2num(handles.source_duration1.String);
        period = str2num(handles.source_period1.String);
        type = handles.source_type1.String{handles.source_type1.Value};
        dt = kgrid.dt;
        t = kgrid.t_array*1e6;
        t0 = find(t >= delay,1);
        t_dur = find(t >= delay+duration,1);
        source.p_mask = zeros(kgrid.Nx,kgrid.Ny);
        for i = 1:length(x_pos)
            source.p_mask(x_pos(i),y_pos(i)) = 1;
        end
        source.p = zeros(1,length(t));
        switch type
            case 'Pulse'
                source.p(t0:t_dur) = magnitude;
            case 'Biphasic'
                t_dur2 = round((t_dur-t0)/2);
                source.p(t0:t_dur) = magnitude;
                source.p(t_dur-t_dur2:t_dur) = -magnitude;
            case 'Sin'
                source.p(t0:t_dur) = magnitude * sin(2*pi*t(1:t_dur-t0+1)/period);
            case 'Cos'
                source.p(t0:t_dur) = magnitude * cos(2*pi*t(1:t_dur-t0+1)/period);
        end
        source_1 = zeros(kgrid.Nx,kgrid.Ny);
    case 'Transducer'
        trans_name = handles.trans_menu.String{handles.trans_menu.Value}
    switch trans_name
        case '1MHz'
            %Get waveform characteristics
            magnitude = str2num(get(handles.source_mag1,'String'));
            delay = str2num(get(handles.source_delay1,'String'));
            duration = str2num(handles.source_duration1.String);
            period = str2num(handles.source_period1.String);
            type = handles.source_type1.String{handles.source_type1.Value};
            dt = kgrid.dt;
            t = kgrid.t_array*1e6;
            t0 = find(t >= delay,1);
            t_dur = find(t >= delay+duration,1);
            source.p = zeros(1,length(t));
            type = handles.source_type1.String{handles.source_type1.Value};
            switch type
                case 'Pulse'
                    source.p(t0:t_dur) = magnitude;
                case 'Biphasic'
                    t_dur2 = round((t_dur-t0)/2);
                    source.p(t0:t_dur) = magnitude;
                    source.p(t_dur-t_dur2:t_dur) = -magnitude;
                case 'Sin'
                    source.p(t0:t_dur) = magnitude * sin(2*pi*t(1:t_dur-t0+1)/period);
                case 'Cos'
                    source.p(t0:t_dur) = magnitude * cos(2*pi*t(1:t_dur-t0+1)/period);
            end
            source.p = filterTimeSeries(kgrid,evalin('base','medium'),source.p);
            %Get transducer characteristics
            arc_pos = str2num(get(handles.trans_arc_pos,'String'));
            radius = str2double(get(handles.trans_radius,'String'));
            diameter = str2double(get(handles.trans_diameter,'String'));
            focus = str2num(get(handles.trans_focus,'String'));
            source.p_mask = makeArc([kgrid.Nx, kgrid.Ny], arc_pos, radius, diameter, focus);
            source_1 = zeros(kgrid.Nx,kgrid.Ny);
        case 'P4-2'
            %get waveform and time grid
             type = handles.source_type1.String{handles.source_type1.Value};
            dt = kgrid.dt;
            t = kgrid.t_array*1e6;
            sampling_freq = 1/kgrid.dt;
            delay = str2num(get(handles.source_delay1,'String'));
            t0 = find(t >= delay,1);
            %get transducer characterstics
            steering_angle = str2double(get(handles.trans_steer,'String'));
            element_spacing = str2double(get(handles.trans_x_kerf,'String'))/1e3;
            burst_freq = str2double(get(handles.trans_freq,'String'))*1e6;
            burst_cycles = str2double(get(handles.trans_cycles,'String'));
            num_elements = str2double(get(handles.trans_elements,'String'));
            x_vec = kgrid.x_vec-min(kgrid.x_vec);
            for i = 1:num_elements
                [~,element_index(i)] = min(abs(x_vec - i*element_spacing));
            end
            
            % use geometric beam forming to calculate the tone burst offsets for each
            % transducer element based on the element index
            medium = evalin('base','medium');
            element_space = kgrid.dx * (-(num_elements - 1)/2:(num_elements - 1)/2);
            burst_offset = t0 + element_space * sin(steering_angle * pi/180)...
                / (medium.sound_speed * kgrid.dt);
            
            % Create tone burst signals
            source.p = toneBurst(sampling_freq, burst_freq, burst_cycles,...
                'SignalOffset',burst_offset);
            source_1 = zeros(kgrid.Nx,kgrid.Ny);
            
            %Create mask
               x0 = str2double(get(handles.trans_x0,'String'));
            element_index = element_index + x0 - round(max(element_index)/2);
            y0 = str2double(get(handles.trans_y0,'String'));
            source.p_mask = zeros(kgrid.Nx,kgrid.Ny);
            source.p_mask(element_index,y0) = 1;
    end
            
end

sourcemenu2 = handles.source_menu2.String{handles.source_menu2.Value};
switch sourcemenu2
    case 'none'
        source_2 = zeros(kgrid.Nx,kgrid.Ny);
    case 'Disk'
        disc_magnitude = str2num(get(handles.source_mag2,'String'));
        disc_x_pos = str2num(get(handles.source_xpos2,'String'));
        disc_y_pos = str2num(get(handles.source_ypos2,'String'));
        disc_rad = str2num(get(handles.source_rad2,'String'));
        source_2 = disc_magnitude * makeDisc(kgrid.Nx, kgrid.Ny, disc_x_pos, disc_y_pos, disc_rad);
end

sourcemenu3 = handles.source_menu3.String{handles.source_menu3.Value};
switch sourcemenu2
    case 'none'
        source_3 = zeros(kgrid.Nx,kgrid.Ny);
    case 'Disk'
        disc_magnitude = str2num(get(handles.source_mag3,'String'));
        disc_x_pos = str2num(get(handles.source_xpos3,'String'));
        disc_y_pos = str2num(get(handles.source_ypos3,'String'));
        disc_rad = str2num(get(handles.source_rad3,'String'));
        source_3 = disc_magnitude * makeDisc(kgrid.Nx, kgrid.Ny, disc_x_pos, disc_y_pos, disc_rad);
end
       
source.p0 = source_1 + source_2 + source_3; %Initial Pressure
assignin('base','source',source);
% hObject    handle to make_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_sensor.
function make_sensor_Callback(hObject, eventdata, handles)
rad = str2num(get(handles.sensor_rad,'String'));
pts = str2num(get(handles.sensor_pts,'String'));

if ismember('kgrid',evalin('base','who'))
    mode_string = get(handles.sensor_mode,'String');
    mode_val = get(handles.sensor_mode,'Value');
    mode = mode_string{mode_val};
    switch mode
        case 'Cart_Circle'
            sensor_rad = rad*1e-3; % [mm]
            sensor_pts = pts; % grid pts
            xpos = str2num(get(handles.sensor_xpos,'String'))*1e-3; %mm from 0,0
            ypos = str2num(get(handles.sensor_ypos,'String'))*1e-3;
            sensor.mask  = makeCartCircle(sensor_rad, sensor_pts,[xpos ypos]);
        case 'Circle'
            xpos = str2num(get(handles.sensor_xpos,'String')); %grid pt
            ypos = str2num(get(handles.sensor_ypos,'String')); %grid pt
            rad = str2num(get(handles.sensor_rad,'String')); %grid pt
            xpts = str2num(get(handles.sensor_xpts,'String')); %grid pt
            ypts = str2num(get(handles.sensor_ypts,'String')); %grid pt
            sensor_arc_angle = 2*pi;
            sensor.mask = makeCircle(xpts,ypts,xpos,ypos,rad,sensor_arc_angle);
        case 'Rect'
            xpos = str2num(get(handles.sensor_xpos,'String')); %grid pt
            ypos = str2num(get(handles.sensor_ypos,'String')); %grid pt
            sensor.mask = [xpos(1), ypos(1), xpos(2), ypos(2)];
        case 'Sphere'
            kgrid = evalin('base','kgrid');
              xpos = str2num(get(handles.sensor_xpos,'String')); %grid pt
            ypos = str2num(get(handles.sensor_ypos,'String')); %grid pt
            zpos = str2num(get(handles.sensor_zpos,'String')); %grid pt
            radius = str2double(handles.sensor_rad.String);
            sensor.mask = makeSphere(kgrid.Nx,kgrid.Ny,kgrid.Nz,radius,1);
            
    end
    rec_p = {'p'};
    if handles.sensor_velocity.Value
        rec_u = {'u'}
    else
        rec_u = {};
    end
    if handles.sensor_final.Value
        if handles.sensor_velocity.Value
        rec_f = {'p_final','u_final'}
        else 
            rec_f = {'p_final'};
        end
    else 
        rec_f = {}; 
    end
    sensor.record = [rec_p, rec_u, rec_f];
    assignin('base','sensor',sensor);
else
    errordlg('Need to create Kgrid before sensor')  
end
% hObject    handle to make_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sensor_rad_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_rad as text
%        str2double(get(hObject,'String')) returns contents of sensor_rad as a double


% --- Executes during object creation, after setting all properties.
function sensor_rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_pts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_pts as text
%        str2double(get(hObject,'String')) returns contents of sensor_pts as a double


% --- Executes during object creation, after setting all properties.
function sensor_pts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sensor_mode.
function sensor_mode_Callback(hObject, eventdata, handles)
switch hObject.String{hObject.Value};
    case 'Cart_Circle'
        set(handles.sensor_rad,'Visible','on');
        set(handles.sensor_rad_text,'Visible','on');
        set(handles.sensor_pts,'Visible','on');
        set(handles.sensor_pts_text,'Visible','on');
        set(handles.sensor_xpos_text,'Visible','on');
        set(handles.sensor_xpos,'Visible','on');
        set(handles.sensor_ypos,'Visible','on');
        set(handles.sensor_ypos_text,'Visible','on');
        set(handles.center_sensor,'Visible','on');
        set(handles.sensor_xpts,'Visible','off');
        set(handles.sensor_xpts_text,'Visible','off');
        set(handles.sensor_ypts,'Visible','off');
        set(handles.sensor_ypts_text,'Visible','off');
        set(handles.sensor_pts_text,'String','Pts');
        set(handles.sensor_rad_text,'String','Radius [mm]');
    case 'Circle'
           set(handles.sensor_rad,'Visible','on');
        set(handles.sensor_rad_text,'Visible','on');
        set(handles.sensor_pts,'Visible','off');
        set(handles.sensor_pts_text,'Visible','off');
        set(handles.sensor_xpos_text,'Visible','on');
        set(handles.sensor_xpos,'Visible','on');    
        set(handles.sensor_ypos,'Visible','on');
        set(handles.sensor_ypos_text,'Visible','on');
        set(handles.center_sensor,'Visible','on');
        set(handles.sensor_xpts,'Visible','on');
        set(handles.sensor_xpts_text,'Visible','on');
        set(handles.sensor_ypts,'Visible','on');
        set(handles.sensor_ypts_text,'Visible','on');
        set(handles.sensor_rad_text,'String','Radius [pts]');
end
% hObject    handle to sensor_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sensor_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sensor_mode


% --- Executes during object creation, after setting all properties.
function sensor_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_xpos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_xpos as text
%        str2double(get(hObject,'String')) returns contents of sensor_xpos as a double


% --- Executes during object creation, after setting all properties.
function sensor_xpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_ypos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_ypos as text
%        str2double(get(hObject,'String')) returns contents of sensor_ypos as a double


% --- Executes during object creation, after setting all properties.
function sensor_ypos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in center_sensor.
function center_sensor_Callback(hObject, eventdata, handles)
if ismember('kgrid',evalin('base','who'))
    switch handles.sensor_mode.String{handles.sensor_mode.Value};
        case 'Cart_Circle'
              set(handles.sensor_xpos,'String',0);
    set(handles.sensor_ypos,'String',0);    
        case 'Circle'
    kgrid = evalin('base','kgrid');
    Ny = kgrid.Ny;
    Nx = kgrid.Nx;
    x = ceil(Nx/2);
    y = ceil(Ny/2);
    set(handles.sensor_xpos,'String',x);
    set(handles.sensor_ypos,'String',y);
    end
else
    errordlg('Need to create Kgrid before sensor')  
end
% hObject    handle to center_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reorder.
function reorder_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sensor = evalin('base','sensor');
sensor_data = evalin('base','sensor_data');
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);
assignin('base','sensor_data',sensor_data_reordered);
% hObject    handle to reorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sensor_ypts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_ypts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_ypts as text
%        str2double(get(hObject,'String')) returns contents of sensor_ypts as a double


% --- Executes during object creation, after setting all properties.
function sensor_ypts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_ypts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_xpts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_xpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_xpts as text
%        str2double(get(hObject,'String')) returns contents of sensor_xpts as a double


% --- Executes during object creation, after setting all properties.
function sensor_xpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_xpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu1.
function source_menu1_Callback(hObject, eventdata, handles)
switch hObject.String{hObject.Value}
    case 'none'
        set(handles.source_text1,'Visible','off')
        set(handles.source_text2,'Visible','off')
        set(handles.source_text3,'Visible','off')
        set(handles.source_text4,'Visible','off')
        set(handles.source_text13,'Visible','off')
        set(handles.source_text14,'Visible','off')
        set(handles.source_text15,'Visible','off')
        set(handles.source_text16,'Visible','off')
        set(handles.source_mag1,'Visible','off')
        set(handles.source_xpos1,'Visible','off')
        set(handles.source_ypos1,'Visible','off')
        set(handles.source_rad1,'Visible','off')
        set(handles.source_delay1,'Visible','off')
        set(handles.source_duration1,'Visible','off')
        set(handles.source_type1,'Visible','off')
        set(handles.source_period1,'Visible','off')
        set(handles.us_transducer,'Visible','off');
    case 'Disk'
        set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','on')
        set(handles.source_text13,'Visible','off')
        set(handles.source_text14,'Visible','off')
        set(handles.source_text15,'Visible','off')
         set(handles.source_text16,'Visible','off')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','on')
        set(handles.source_delay1,'Visible','off')
        set(handles.source_duration1,'Visible','off')
        set(handles.source_type1,'Visible','off')
        set(handles.source_period1,'Visible','off')
        set(handles.us_transducer,'Visible','off');
    case 'Time-Varying Point'
        set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','off')
        set(handles.source_text13,'Visible','on')
        set(handles.source_text14,'Visible','on')
        set(handles.source_text15,'Visible','on')
        set(handles.source_text16,'Visible','on')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','off')
        set(handles.source_delay1,'Visible','on')
        set(handles.source_duration1,'Visible','on')
        set(handles.source_type1,'Visible','on')
        set(handles.source_period1,'Visible','on')
        set(handles.us_transducer,'Visible','off');
    case 'Transducer'
          set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','on')
        set(handles.source_text13,'Visible','on')
        set(handles.source_text14,'Visible','on')
        set(handles.source_text15,'Visible','on')
        set(handles.source_text16,'Visible','on')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','on')
        set(handles.source_delay1,'Visible','on')
        set(handles.source_duration1,'Visible','on')
        set(handles.source_type1,'Visible','on')
        set(handles.source_period1,'Visible','on')
        set(handles.source_menu2,'Value',1);
        set(handles.source_menu3,'Value',1);
        set(handles.us_transducer,'Visible','on');
        
end
% hObject    handle to source_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu1


% --- Executes during object creation, after setting all properties.
function source_menu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu2.
function source_menu2_Callback(hObject, eventdata, handles)
switch hObject.String{hObject.Value}
    case 'none'
        set(handles.source_text5,'Visible','off')
        set(handles.source_text6,'Visible','off')
        set(handles.source_text7,'Visible','off')
        set(handles.source_text8,'Visible','off')
        set(handles.source_mag2,'Visible','off')
        set(handles.source_xpos2,'Visible','off')
        set(handles.source_ypos2,'Visible','off')
        set(handles.source_rad2,'Visible','off')
    case 'Disk'
          set(handles.source_text5,'Visible','on')
        set(handles.source_text6,'Visible','on')
        set(handles.source_text7,'Visible','on')
        set(handles.source_text8,'Visible','on')
        set(handles.source_mag2,'Visible','on')
        set(handles.source_xpos2,'Visible','on')
        set(handles.source_ypos2,'Visible','on')
        set(handles.source_rad2,'Visible','on')
end
% hObject    handle to source_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu2


% --- Executes during object creation, after setting all properties.
function source_menu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu3.
function source_menu3_Callback(hObject, eventdata, handles)
switch hObject.String{hObject.Value}
    case 'none'
        set(handles.source_text9,'Visible','off')
        set(handles.source_text10,'Visible','off')
        set(handles.source_text11,'Visible','off')
        set(handles.source_text12,'Visible','off')
        set(handles.source_mag3,'Visible','off')
        set(handles.source_xpos3,'Visible','off')
        set(handles.source_ypos3,'Visible','off')
        set(handles.source_rad3,'Visible','off')
    case 'Disk'
          set(handles.source_text9,'Visible','on')
        set(handles.source_text10,'Visible','on')
        set(handles.source_text11,'Visible','on')
        set(handles.source_text12,'Visible','on')
        set(handles.source_mag3,'Visible','on')
        set(handles.source_xpos3,'Visible','on')
        set(handles.source_ypos3,'Visible','on')
        set(handles.source_rad3,'Visible','on')
end
% hObject    handle to source_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu3


% --- Executes during object creation, after setting all properties.
function source_menu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag1_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag1 as text
%        str2double(get(hObject,'String')) returns contents of source_mag1 as a double


% --- Executes during object creation, after setting all properties.
function source_mag1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos1_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos1 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos1 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos1_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos1 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos1 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad1_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad1 as text
%        str2double(get(hObject,'String')) returns contents of source_rad1 as a double


% --- Executes during object creation, after setting all properties.
function source_rad1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag2_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag2 as text
%        str2double(get(hObject,'String')) returns contents of source_mag2 as a double


% --- Executes during object creation, after setting all properties.
function source_mag2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos2_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos2 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos2 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos2_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos2 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos2 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad2_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad2 as text
%        str2double(get(hObject,'String')) returns contents of source_rad2 as a double


% --- Executes during object creation, after setting all properties.
function source_rad2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag3_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag3 as text
%        str2double(get(hObject,'String')) returns contents of source_mag3 as a double


% --- Executes during object creation, after setting all properties.
function source_mag3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos3_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos3 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos3 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos3_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos3 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos3 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad3_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad3 as text
%        str2double(get(hObject,'String')) returns contents of source_rad3 as a double


% --- Executes during object creation, after setting all properties.
function source_rad3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sensor_rad_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_rad_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in record_movie.
function record_movie_Callback(hObject, eventdata, handles)
% hObject    handle to record_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of record_movie



function sim_scale_Callback(hObject, eventdata, handles)
% hObject    handle to sim_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_scale as text
%        str2double(get(hObject,'String')) returns contents of sim_scale as a double


% --- Executes during object creation, after setting all properties.
function sim_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sim_freq_Callback(hObject, eventdata, handles)
% hObject    handle to sim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_freq as text
%        str2double(get(hObject,'String')) returns contents of sim_freq as a double


% --- Executes during object creation, after setting all properties.
function sim_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sim_fps_Callback(hObject, eventdata, handles)
% hObject    handle to sim_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_fps as text
%        str2double(get(hObject,'String')) returns contents of sim_fps as a double


% --- Executes during object creation, after setting all properties.
function sim_fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sim_mesh.
function sim_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to sim_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_mesh


% --- Executes on button press in sim_mask.
function sim_mask_Callback(hObject, eventdata, handles)
% hObject    handle to sim_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_mask


% --- Executes on button press in sim_pml.
function sim_pml_Callback(hObject, eventdata, handles)
% hObject    handle to sim_pml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_pml


% --- Executes on button press in sensor_velocity.
function sensor_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_velocity


% --- Executes on selection change in medium_menu.
function medium_menu_Callback(hObject, eventdata, handles)
V = hObject.Value
switch V
    case 1
        set(handles.medium_layers,'Visible','off')
    case 2
        set(handles.medium_layers,'Visible','on');
end
% hObject    handle to medium_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns medium_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from medium_menu


% --- Executes during object creation, after setting all properties.
function medium_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_ss_Callback(hObject, eventdata, handles)
% hObject    handle to medium_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_ss as text
%        str2double(get(hObject,'String')) returns contents of medium_ss as a double


% --- Executes during object creation, after setting all properties.
function medium_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_density_Callback(hObject, eventdata, handles)
% hObject    handle to medium_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_density as text
%        str2double(get(hObject,'String')) returns contents of medium_density as a double


% --- Executes during object creation, after setting all properties.
function medium_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_alpha_coeff_Callback(hObject, eventdata, handles)
% hObject    handle to medium_alpha_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_alpha_coeff as text
%        str2double(get(hObject,'String')) returns contents of medium_alpha_coeff as a double


% --- Executes during object creation, after setting all properties.
function medium_alpha_coeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_alpha_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_alpha_power_Callback(hObject, eventdata, handles)
% hObject    handle to medium_alpha_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_alpha_power as text
%        str2double(get(hObject,'String')) returns contents of medium_alpha_power as a double


% --- Executes during object creation, after setting all properties.
function medium_alpha_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_alpha_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data_menu.
function data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_menu


% --- Executes during object creation, after setting all properties.
function data_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ly1_rnd.
function ly1_rnd_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_rnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly1_rnd



function ly3_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_yN as text
%        str2double(get(hObject,'String')) returns contents of ly3_yN as a double


% --- Executes during object creation, after setting all properties.
function ly3_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_yN as text
%        str2double(get(hObject,'String')) returns contents of ly2_yN as a double


% --- Executes during object creation, after setting all properties.
function ly2_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_yN as text
%        str2double(get(hObject,'String')) returns contents of ly1_yN as a double


% --- Executes during object creation, after setting all properties.
function ly1_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly3_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly3_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_xN as text
%        str2double(get(hObject,'String')) returns contents of ly3_xN as a double


% --- Executes during object creation, after setting all properties.
function ly3_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly3_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly3_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly2_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly2_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_xN as text
%        str2double(get(hObject,'String')) returns contents of ly2_xN as a double


% --- Executes during object creation, after setting all properties.
function ly2_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly2_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly2_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly1_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly1_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_xN as text
%        str2double(get(hObject,'String')) returns contents of ly1_xN as a double


% --- Executes during object creation, after setting all properties.
function ly1_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly1_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly1_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ly1_on.
function ly1_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly1_on


% --- Executes on button press in ly2_on.
function ly2_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly2_on


% --- Executes on button press in ly3_on.
function ly3_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly3_on



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_ac as text
%        str2double(get(hObject,'String')) returns contents of ly3_ac as a double


% --- Executes during object creation, after setting all properties.
function ly3_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_den as text
%        str2double(get(hObject,'String')) returns contents of ly3_den as a double


% --- Executes during object creation, after setting all properties.
function ly3_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_ss as text
%        str2double(get(hObject,'String')) returns contents of ly3_ss as a double


% --- Executes during object creation, after setting all properties.
function ly3_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_ac as text
%        str2double(get(hObject,'String')) returns contents of ly2_ac as a double


% --- Executes during object creation, after setting all properties.
function ly2_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_den as text
%        str2double(get(hObject,'String')) returns contents of ly2_den as a double


% --- Executes during object creation, after setting all properties.
function ly2_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_ss as text
%        str2double(get(hObject,'String')) returns contents of ly2_ss as a double


% --- Executes during object creation, after setting all properties.
function ly2_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_ac as text
%        str2double(get(hObject,'String')) returns contents of ly1_ac as a double


% --- Executes during object creation, after setting all properties.
function ly1_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_den as text
%        str2double(get(hObject,'String')) returns contents of ly1_den as a double


% --- Executes during object creation, after setting all properties.
function ly1_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_ss as text
%        str2double(get(hObject,'String')) returns contents of ly1_ss as a double


% --- Executes during object creation, after setting all properties.
function ly1_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_medium.
function plot_medium_Callback(hObject, eventdata, handles)
medium = evalin('base','medium');
axes(handles.axes2)
imagesc(medium.sound_speed);
xlabel('Y Position');
ylabel('X Position');
title('Sound Speed');
colormap('jet');
colorbar;

axes(handles.axes3)
hold off;
imagesc(medium.density);
xlabel('Y Position');
ylabel('X Position');
title('Density');
colormap('jet');
colorbar;
% hObject    handle to plot_medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function source_delay1_Callback(hObject, eventdata, handles)
% hObject    handle to source_delay1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_delay1 as text
%        str2double(get(hObject,'String')) returns contents of source_delay1 as a double


% --- Executes during object creation, after setting all properties.
function source_delay1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_delay1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_duration1_Callback(hObject, eventdata, handles)
% hObject    handle to source_duration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_duration1 as text
%        str2double(get(hObject,'String')) returns contents of source_duration1 as a double


% --- Executes during object creation, after setting all properties.
function source_duration1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_duration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_type1.
function source_type1_Callback(hObject, eventdata, handles)
% hObject    handle to source_type1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_type1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_type1


% --- Executes during object creation, after setting all properties.
function source_type1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_type1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_period1_Callback(hObject, eventdata, handles)
% hObject    handle to source_period1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_period1 as text
%        str2double(get(hObject,'String')) returns contents of source_period1 as a double


% --- Executes during object creation, after setting all properties.
function source_period1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_period1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sensor_final.
function sensor_final_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_final



function display_clims_Callback(hObject, eventdata, handles)
% hObject    handle to display_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_clims as text
%        str2double(get(hObject,'String')) returns contents of display_clims as a double


% --- Executes during object creation, after setting all properties.
function display_clims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trans_menu.
function trans_menu_Callback(hObject, eventdata, handles)
trans = handles.trans_menu.String{handles.trans_menu.Value};
kgrid = evalin('base','kgrid');
switch trans
    case '1MHz'
        set(handles.source_period1,'String',1)
        set(handles.source_duration1,'String',2)
        set(handles.source_delay1,'String',0)
        x = kgrid.Nx/2;
        y = kgrid.Ny-4;
        set(handles.trans_arc_pos,'String',num2str([x y]))
        r = 68/kgrid.dx/1e3;
        d = round((30/kgrid.dx)/1e3)+1;
        set(handles.trans_radius,'String',r);
        set(handles.trans_diameter,'String',d);
        f = [x round(68/kgrid.dx/1e3)-4];
        set(handles.trans_focus,'String',num2str(f));
    case 'P4-1'
    case 'P4-2'
        x = kgrid.Nx/2;
        set(handles.trans_x0,'String',x); %grid
        y = kgrid.Ny-4;
        set(handles.trans_y0,'String',y); %grid
        xkerf = 0.3/kgrid.dx/1e3; 
        length = 0.4831/kgrid.dx/1e3;
        set(handles.trans_x_kerf,'String',num2str(xkerf)); %mm
        set(handles.trans_width,'String',num2str(length));
        set(handles.trans_length,'String',num2str(length));
        set(handles.trans_cycles,'String',2);
        set(handles.trans_freq,'String',3); %MHz
        set(handles.trans_steer,'String',0); %degrees
        set(handles.trans_el_focus,'String',60); %mm
        set(handles.trans_lens_focus,'String',60) %mm
        set(handles.trans_active_elements,'String','1:64');
    case 'H235'
    case 'h247'
end
% hObject    handle to trans_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trans_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trans_menu


% --- Executes during object creation, after setting all properties.
function trans_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_x0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_x0 as text
%        str2double(get(hObject,'String')) returns contents of trans_x0 as a double


% --- Executes during object creation, after setting all properties.
function trans_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_y0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_y0 as text
%        str2double(get(hObject,'String')) returns contents of trans_y0 as a double


% --- Executes during object creation, after setting all properties.
function trans_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_arc_pos_Callback(hObject, eventdata, handles)
% hObject    handle to trans_arc_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_arc_pos as text
%        str2double(get(hObject,'String')) returns contents of trans_arc_pos as a double


% --- Executes during object creation, after setting all properties.
function trans_arc_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_arc_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_radius_Callback(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_radius as text
%        str2double(get(hObject,'String')) returns contents of trans_radius as a double


% --- Executes during object creation, after setting all properties.
function trans_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to trans_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_diameter as text
%        str2double(get(hObject,'String')) returns contents of trans_diameter as a double


% --- Executes during object creation, after setting all properties.
function trans_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function record_filename_Callback(hObject, eventdata, handles)
% hObject    handle to record_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of record_filename as text
%        str2double(get(hObject,'String')) returns contents of record_filename as a double


% --- Executes during object creation, after setting all properties.
function record_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function record_dir_Callback(hObject, eventdata, handles)
% hObject    handle to record_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of record_dir as text
%        str2double(get(hObject,'String')) returns contents of record_dir as a double


% --- Executes during object creation, after setting all properties.
function record_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_elements_Callback(hObject, eventdata, handles)
% hObject    handle to trans_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_elements as text
%        str2double(get(hObject,'String')) returns contents of trans_elements as a double


% --- Executes during object creation, after setting all properties.
function trans_elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_x_kerf_Callback(hObject, eventdata, handles)
% hObject    handle to trans_x_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_x_kerf as text
%        str2double(get(hObject,'String')) returns contents of trans_x_kerf as a double


% --- Executes during object creation, after setting all properties.
function trans_x_kerf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_x_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_steer_Callback(hObject, eventdata, handles)
% hObject    handle to trans_steer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_steer as text
%        str2double(get(hObject,'String')) returns contents of trans_steer as a double


% --- Executes during object creation, after setting all properties.
function trans_steer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_steer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_freq_Callback(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_freq as text
%        str2double(get(hObject,'String')) returns contents of trans_freq as a double


% --- Executes during object creation, after setting all properties.
function trans_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_cycles_Callback(hObject, eventdata, handles)
% hObject    handle to trans_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_cycles as text
%        str2double(get(hObject,'String')) returns contents of trans_cycles as a double


% --- Executes during object creation, after setting all properties.
function trans_cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sim_3D.
function sim_3D_Callback(hObject, eventdata, handles)
% hObject    handle to sim_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_3D


% --- Executes on button press in make_trans.
function make_trans_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
medium = evalin('base','medium');
tw = evalin('base','tw');
transducer.number_elements = str2double(handles.trans_elements.String);    % total number of transducer elements
e_width = str2double(handles.trans_width.String); 
e_width = round(e_width/kgrid.dy/1e3)
transducer.element_width = e_width;       % width of each element [grid points]
e_length = str2double(handles.trans_length.String);
e_length = round(e_length/kgrid.dz/1e3);
transducer.element_length = e_length;    % length of each element [grid points]
e_kerf = str2double(handles.trans_x_kerf.String);
e_kerf = round(e_kerf/kgrid.dy/1e3);
transducer.element_spacing = e_kerf;    % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;           % radius of curvature of the transducer [m]
% Denotes center of transducer position
x0 = str2double(handles.trans_x0.String);
y0 = str2double(handles.trans_y0.String);
z0 = str2double(handles.trans_z0.String);
transducer.position = round([x0,y0,z0]);

%Beamforming delay properities
%transducer.sound_speed = medium.sound_speed; %[m/s]
transducer.sound_speed = medium.sound_speed(1,1); %[m/s]
% focus_loc = str2num(handles.trans_focus.String);
% focus = [kgrid.x_vec(focus_loc(1)) kgrid.y_vec(focus_loc(2)) kgrid.z_vec(focus_loc(3))];
transducer.focus_distance = str2double(handles.trans_lens_focus.String); %[m]
transducer.elevation_focus_distance = str2double(handles.trans_el_focus.String);
transducer.steering_angle = str2double(handles.trans_steer.String);

% Apodization and active elements
transducer.transmit_apodization = handles.t_apod.String{handles.t_apod.Value};
transducer.receive_apodization = handles.r_apod.String{handles.r_apod.Value};
transducer.active_elements = zeros(transducer.number_elements,1);
active = str2num(handles.trans_active_elements.String);
transducer.active_elements(active) = 1;

% Create transducer
tw = evalin('base','tw');
transducer.input_signal = tw;
transducer = kWaveTransducer(kgrid, transducer);
assignin('base','transducer',transducer);


% hObject    handle to make_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_trans.
function show_trans_Callback(hObject, eventdata, handles)
% hObject    handle to show_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_tw.
function show_tw_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
medium = evalin('base','medium');
strength = str2double(handles.source_mag1.String)*1e6;
freq = str2double(handles.trans_freq.String)*1e6;
cycles = str2double(handles.trans_cycles.String);
signal = toneBurst(1/kgrid.dt,freq,cycles);
signal = signal .*strength;
axes(handles.axes2)
plot(kgrid.t_array(1:length(signal)),signal);
xlabel('Time [us]')
ylabel('Particle Velocity [m/s]')
axes(handles.axes3)
dt = kgrid.dt;
maxf = 1/dt/2;
f_axis = linspace(0,maxf,length(signal)/2);
Signal = fft(signal);
plot(f_axis/1e6,abs(Signal(1:length(f_axis))))
xlabel('Frequency [MHz]')
ylabel('Amplitude [au]')
assignin('base','tw',signal);

% hObject    handle to show_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function trans_width_Callback(hObject, eventdata, handles)
% hObject    handle to trans_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_width as text
%        str2double(get(hObject,'String')) returns contents of trans_width as a double


% --- Executes during object creation, after setting all properties.
function trans_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_length_Callback(hObject, eventdata, handles)
% hObject    handle to trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_length as text
%        str2double(get(hObject,'String')) returns contents of trans_length as a double


% --- Executes during object creation, after setting all properties.
function trans_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wave_type.
function wave_type_Callback(hObject, eventdata, handles)
% hObject    handle to wave_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wave_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wave_type


% --- Executes during object creation, after setting all properties.
function wave_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in center_trans.
function center_trans_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid')
y = kgrid.Ny/2;
x = kgrid.Nx/2;
set(handles.trans_x0,'String',x)
set(handles.trans_y0,'String',y)
if kgrid.Nz > 0
    z = kgrid.Nz/2;
    set(handles.trans_z0,'String',z);
end
% hObject    handle to center_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function trans_z0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_z0 as text
%        str2double(get(hObject,'String')) returns contents of trans_z0 as a double


% --- Executes during object creation, after setting all properties.
function trans_z0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_lens_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_lens_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_lens_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_lens_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_lens_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_lens_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_el_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_el_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_el_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_el_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_el_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_el_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in t_apod.
function t_apod_Callback(hObject, eventdata, handles)
% hObject    handle to t_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns t_apod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t_apod


% --- Executes during object creation, after setting all properties.
function t_apod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in r_apod.
function r_apod_Callback(hObject, eventdata, handles)
% hObject    handle to r_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns r_apod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from r_apod


% --- Executes during object creation, after setting all properties.
function r_apod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_active_elements_Callback(hObject, eventdata, handles)
% hObject    handle to trans_active_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_active_elements as text
%        str2double(get(hObject,'String')) returns contents of trans_active_elements as a double


% --- Executes during object creation, after setting all properties.
function trans_active_elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_active_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_zpos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_zpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_zpos as text
%        str2double(get(hObject,'String')) returns contents of sensor_zpos as a double


% --- Executes during object creation, after setting all properties.
function sensor_zpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_zpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
