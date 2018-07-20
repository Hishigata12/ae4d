function varargout = Stitch(varargin)
% Macro to stitch or average two images together into a single image.
% Brings up a GUI to allow use to select two images from listboxes,
% select whether they want them stitched side-by-side, or over-and-under,
% where to place smaller images in the event of size mismatches, and allows
% them to rotate one or both images by 90 or 180 degrees.  It can process a
% single pair of images, or batch process a number of selected pairs.
%
% Author: Mark Hayworth, The Procter & Gamble Company, Cincinnati Ohio
% Date:   May 2007 - October 2009
%
% You are permitted to make changes to this macro and redistribute it as
% long as you keep my name as author and add your on the "Revised by" line.
%
% Revised by:
%
% Stitch M-file for Stitch.fig
%      Stitch by itself, creates a new Stitch or raises the existing
%      singleton*.
%
%      H = Stitch returns the handle to a new Stitch or the handle to
%      the existing singleton*.
%
%      Stitch('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Stitch.M with the given input arguments.
%
%      Stitch('Property','Value',...) creates a new Stitch or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Stitch_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Stitch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Stitch

% Last Modified by GUIDE v2.5 08-Nov-2007 13:58:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Stitch_OpeningFcn, ...
                   'gui_OutputFcn',  @Stitch_OutputFcn, ...
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

% --- Executes just before Stitch is made visible.
function Stitch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Stitch (see VARARGIN)

% Choose default command line output for Stitch
handles.output = hObject;

%=====================================================================
% --- My Startup Code --------------------------------------------------
    % Initialize image folder that the listbox is looking at.
    % Clear old stuff from console.
	clc;
	
	% Change the current folder to the folder of this m-file.
	cd(fileparts(which(mfilename)));
	
	handles.macroFolder = cd;

	% Load up the initial values from the mat file.
	handles = RecallParameters(handles);
	bdfunction = get(handles.axesLeftImage, 'ButtonDownFcn');
    
    % Load splash images.
	LoadSplashImages(handles);
	
	% Enlarge and center window.
	% calculate position in normalized units
	widthFraction = 0.9;
	heightFraction = 0.85;
	positionProperty = [(1 - widthFraction)/2, (1 - heightFraction)/2, widthFraction, heightFraction]; % [left, bottom, width, height]

	% display figure
	set(gcf, 'Units','Normalized', 'Position', positionProperty);



	% !!!! QUIRK workaround.!!!!
	% Make it so that if they click in the image axes, it will execute the
	% button down callback.  Now it won't  -- unless you do this quirk
	% workaround.  For further explanation, see
	% http://www.mathworks.com/support/bugreports/details.html?rp=296908
	%set(axesChildHandle, 'HitTest', 'off');
	% Now double click on the dialog box's background to bring up the
	% property inspector.  Change the WindowButtonDownFcn property so that
	% it says Stitch('Viewer_ButtonDownFcn',gcbo,[],guidata(gcbo))
	% Then do this:
	%set(axesChildHandle, 'ButtonDownFcn', 'Viewer_ButtonDownFcn');
	% !!!! End QUIRK workaround.!!!!

    % Update handles structure
    guidata(hObject, handles);

% --- End of My Startup Code --------------------------------------------------
%=====================================================================


%=====================================================================
% --- Outputs from this function are returned to the command line.
function varargout = Stitch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%=====================================================================
function LoadSplashImages(handles)
	% Load left splash image into left image buffer.
    strSplashFullFileName = [handles.macroFolder '\Stitch1.jpg'];
	leftSplashImage = DisplayImage(handles, 'Left', strSplashFullFileName);
	% Load left splash image into left image buffer.
    strSplashFullFileName = [handles.macroFolder '\Stitch2.jpg'];
	rightSplashImage = DisplayImage(handles, 'Right', strSplashFullFileName);
	% Create the stitched image.
	stitchedImage = [leftSplashImage rightSplashImage];
	% Put stitched image into axesStitchedImage.
	axes(handles.axesStitchedImage);
	imshow(stitchedImage);
	return;	% LoadSplashImages
	
%=====================================================================
% --- Executes on clicking in lstLeftImageList listbox.
% Display image from disk and plots histogram
function lstLeftImageList_Callback(hObject, eventdata, handles)
% hObject    handle to lstLeftImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns lstLeftImageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstLeftImageList

	% Make the image array global so we can bring it up in imtool via the
	% windowbuttondownfcn function if they click on the image in the GUI.
	global leftImage;
	
	axes(handles.axesLeftImage);
	% Change mouse pointer (cursor) to an hourglass.  
	% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.

	% Update the number of images in the Analyze button caption.
	handles = UpdateAnalyzeButtonCaption(handles);

	% Get image name
    Selected = get(handles.lstLeftImageList, 'value');
    % If more than one is selected, bail out.
    if length(Selected) > 1 
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow')
		drawnow;	% Cursor won't change right away unless you do this.
        return;
    end
    % If only one is selected, display it.
    ListOfImageNames = get(handles.lstLeftImageList, 'string');
    baseImageFileName = strcat(cell2mat(ListOfImageNames(Selected)));
    fullImageFileName = [handles.leftImageFolder '\' baseImageFileName];	% Prepend folder.
	set(handles.txtLeftImageName, 'String', baseImageFileName);
	
	[folder, baseImageFileName, extension] = fileparts(fullImageFileName);
	switch lower(extension)
	case '.wav'
		try
			[y, Fs, nbits] = wavread(fullImageFileName);
			plot(y);
			wavplay(y, Fs, 'async');
		catch
			strError = lasterror;
			strErrorMessage = sprintf('MATLAB does not support this type of wave file.\n\n%s', strError.message);
			msgboxw(strErrorMessage);
		end
	case {'.mov', '.wmv', '.asf', '.avi'}
		msgboxw('Video files are not supported by Stitch.');
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		return;
	otherwise
		% Display the image.
		leftImage = DisplayImage(handles, 'Left', fullImageFileName);

		% Display size parameters of the image below the listbox.
		imageSize = size(leftImage);
		strLabel = sprintf('%d x %d', imageSize(1), imageSize(2));
		if length(imageSize) >= 3
			strLabel = [strLabel ' RGB'];
		else
			strLabel = [strLabel ' Monochrome'];
		end
		set(handles.txtLeftImageSize, 'String', strLabel);

		% Rotate it if the checkbox is checked.
		showRotated = get(handles.chkShowRotated, 'Value');
		if showRotated == 1
			% Rotate the right image because they checked the box.
			selectedRotationIndex = get(handles.popLeftRotation,'Value'); % returns selected item from popRightRotation
			switch(selectedRotationIndex)
				case 2
					% 90 clockwise
					leftImage = imrotate(leftImage, -90);
				case 3
					% 90 counter clockwise
					leftImage = imrotate(leftImage, 90);
				case 4
					% 180
					leftImage = imrotate(leftImage, 180);
			end
			% Display the original image rotated, as long as it's not 0 degrees.
			if selectedRotationIndex ~= 1
				axes(handles.axesLeftImage);
				imshow(leftImage);
			end
		end

	end
	     
	% Change mouse pointer (cursor) to an arrow.
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
    guidata(hObject, handles);
    return; % lstImageList_Callback

% --- Executes on selection change in lstRightImageList.
function lstRightImageList_Callback(hObject, eventdata, handles)
% hObject    handle to lstRightImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns lstRightImageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstRightImageList
	axes(handles.axesRightImage);

	% Make the image array global so we can bring it up in imtool via the
	% windowbuttondownfcn function if they click on the image in the GUI.
	global rightImage;

	% Change mouse pointer (cursor) to an hourglass.  
	% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.

	% Update the number of images in the Analyze button caption.
	handles = UpdateAnalyzeButtonCaption(handles);

	% Get image name
    Selected = get(handles.lstRightImageList, 'value');
    % If more than one is selected, bail out.
    if length(Selected) > 1 
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow')
		drawnow;	% Cursor won't change right away unless you do this.
        return;
    end
    % If only one is selected, display it.
    ListOfImageNames = get(handles.lstRightImageList, 'string');
    baseImageFileName = strcat(cell2mat(ListOfImageNames(Selected)));
    fullImageFileName = [handles.rightImageFolder '\' baseImageFileName];	% Prepend folder.
	set(handles.txtRightImageName, 'String', baseImageFileName);
	
	[folder, baseFileName, extension] = fileparts(fullImageFileName);
	switch lower(extension)
	case '.wav'
		try
			[y, Fs, nbits] = wavread(fullImageFileName);
			plot(y);
			wavplay(y, Fs, 'async');
		catch
			strError = lasterror;
			strErrorMessage = sprintf('MATLAB does not support this type of wave file.\n\n%s', strError.message);
			msgboxw(strErrorMessage);
		end
	case {'.mov', '.wmv', '.asf', '.avi'}
		msgboxw('Video files are not supported by Stitch.');
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		return;
	otherwise
		% Display the image.
		rightImage = DisplayImage(handles, 'Right', fullImageFileName);

		% Display size parameters of the image below the listbox.
		imageSize = size(rightImage);
		strLabel = sprintf('%d x %d', imageSize(1), imageSize(2));
		if length(imageSize) >= 3
			strLabel = [strLabel ' RGB'];
		else
			strLabel = [strLabel ' Monochrome'];
		end
		set(handles.txtRightImageSize, 'String', strLabel);

		% Rotate it if the checkbox is checked.
		showRotated = get(handles.chkShowRotated, 'Value');
		if showRotated == 1
			% Rotate the right image if necessary.
			selectedRotationIndex = get(handles.popRightRotation,'Value'); % returns selected item from popRightRotation
			switch(selectedRotationIndex)
				case 2
					% 90 clockwise
					rightImage = imrotate(rightImage, -90);
				case 3
					% 90 counter clockwise
					rightImage = imrotate(rightImage, 90);
				case 4
					% 180
					rightImage = imrotate(rightImage, 180);
			end
			% Display the original image rotated, as long as it's not 0 degrees.
			if selectedRotationIndex ~= 1
				axes(handles.axesRightImage);
				imshow(rightImage);
			end
		end
	end
	     
	% Change mouse pointer (cursor) to an arrow.
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
    guidata(hObject, handles);
    return; % lstImageList_Callback


%=====================================================================
% Reads FullImageFileName from disk into the axesToUse axes.
function imageArray = DisplayImage(handles, strAxesToUse, FullImageFileName)
	if strcmpi(strAxesToUse, 'Left')
		axesToUse = handles.axesLeftImage;
	elseif strcmpi(strAxesToUse, 'Right')
		axesToUse = handles.axesRightImage;
	elseif strcmpi(strAxesToUse, 'Stitched')
		axesToUse = handles.axesStitchedImage;
		set(handles.axesStitchedImage, 'visible', 'on');
	end
	
	% Read from disk into an array.
	imageArray = imread(FullImageFileName);

    % Display image array in a window on the user interface.
    axes(axesToUse);
    handleToImage = imshow(imageArray, []);

	% Set the button down function for the image to be the one we want it
	% to be (the one that will bring up the proper image in imtool).
	if strcmpi(strAxesToUse, 'Left')
		set(handleToImage, 'ButtonDownFcn', @Left_Image_Click_Event);
	elseif strcmpi(strAxesToUse, 'Right')
		set(handleToImage, 'ButtonDownFcn', @Right_Image_Click_Event);
	elseif strcmpi(strAxesToUse, 'Stitched')
		set(handleToImage, 'ButtonDownFcn', @Stitched_Image_Click_Event);
	end

	return; % DisplayImage


%=====================================================================
function handles = UpdateAnalyzeButtonCaption(handles)
	% Find out if they specified Stitch or Blend.
	selectedIndex = get(handles.popMode, 'Value');
	%mode = contents{selectedIndex}
	if selectedIndex == 1
		modeString = 'Stitch';
	else
		modeString = 'Blend';
	end

	% Get what indices were selected on both sides.
    leftSelected = get(handles.lstLeftImageList, 'value');
    rightSelected = get(handles.lstRightImageList, 'value');
	% Find out the number of indices that were selected.
	numberSelectedOnLeft = length(leftSelected);
	numberSelectedOnRight = length(rightSelected);
	if numberSelectedOnLeft == numberSelectedOnRight && numberSelectedOnLeft >= 1
		% The same number has been selected on each side, 
		% and that number is greater than 0.  Enable Stitch button.
		buttonCaption = sprintf('Step 5:  %s ', modeString);
        buttonCaption = {['Step 5:  ' modeString ' ']};   % MATLAB quirk - needs to be cell array to keep trailing spaces.
        buttonCaption = strcat(buttonCaption, num2str(numberSelectedOnLeft));
		if numberSelectedOnLeft == 1 
		    buttonCaption = strcat(buttonCaption, ' image pair');
		else
		    buttonCaption = strcat(buttonCaption, ' image pairs');
		end
        set(handles.btnStitch, 'string', buttonCaption);
        set(handles.btnStitch, 'Enable', 'on');
	else
		% Uneven number of images selected, or no images selected.
		% Disable the Stitch button.
        set(handles.btnStitch, 'string', 'Step 5:  Stitch no images');
        set(handles.btnStitch, 'Enable', 'off');
	end
	return;
		
%=====================================================================
% --- Executes during object creation, after setting all properties.
function lstLeftImageList_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lstLeftImageList (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    return

%=====================================================================
% --- Executes on clicking btnSelectFolder button.
% Asks user to select a directory and then loads up the listbox (via a call
% to LoadImageList)
function btnSelectFolder_Callback(hObject, eventdata, handles)
    % hObject    handle to btnSelectFolder (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    %msgbox(handles.ImageFolder);
	folderToGet = get(handles.radLeftImageFolder, 'Value');
	if folderToGet == 1	% This radio button is selected
		% Specify the left folder.
		returnValue = uigetdir(handles.leftImageFolder,'Select left folder');
		% returnValue will be 0 (a double) if they click cancel.
		% returnValue will be the path (a string) if they clicked OK.
		if returnValue ~= 0
			% Assign the value if they didn't click cancel.
			handles.leftImageFolder = returnValue;
		end
		% Load up the image listbox.
		handles = LoadImageList(handles, 1, 0, 0);
	else % The left radio button is not selected.
		% Specify the right folder, or the stitched folder.
		folderToGet = get(handles.radRightImageFolder, 'Value');
		if folderToGet == 1	% The right radio button is selected
			returnValue = uigetdir(handles.rightImageFolder,'Select right folder');
			% returnValue will be 0 (a double) if they click cancel.
			% returnValue will be the path (a string) if they clicked OK.
			if returnValue ~= 0
				% Assign the value if they didn't click cancel.
				handles.rightImageFolder = returnValue;
			end
			% Load up the image listbox.
			handles = LoadImageList(handles, 0, 1, 0);
		else
			% The Stitched button is selected.
			returnValue = uigetdir(handles.stitchedImageFolder,'Select folder to contain stitched images.');
			% returnValue will be 0 (a double) if they click cancel.
			% returnValue will be the path (a string) if they clicked OK.
			if returnValue ~= 0
				% Assign the value if they didn't click cancel.
				handles.stitchedImageFolder = returnValue;
			end
			% Load up the image listbox.
			handles = LoadImageList(handles, 0, 0, 1);
		end
	end

	% Assign to the labels over the listboxes.
	set(handles.txtLeftFolder, 'string' ,handles.leftImageFolder);
	set(handles.txtRightFolder, 'string' ,handles.rightImageFolder);

	% Save the image folders in our initialization mat file.
	SaveParameters(handles);
	guidata(hObject, handles);
    return;

	
%=====================================================================
% Saves the folders and other user option settings to the mat file.
function handles = SaveParameters(handles)
	% Save the image folders in our ini file.
	lastUsedLeftImageFolder = handles.leftImageFolder;
	lastUsedRightImageFolder = handles.rightImageFolder;
	lastUsedStitchedImageFolder = handles.stitchedImageFolder;
	
	% Get the values of the optional parameters so we can save those too.
	lastUsedOrientationOption = get(handles.popOrientation, 'value');
	lastUsedSizeMismatchOption = get(handles.popSizeMismatch, 'value');
	lastUsedPauseOption = get(handles.chkPauseAfterImage, 'value');
	lastUsedOverwriteOption = get(handles.chkPromptToOverwrite, 'value');
	lastUsedSaveFileNameOption = get(handles.chkPromptForFileName, 'value');
	
	% Write everything out.
	save('Stitch.mat', 'lastUsedLeftImageFolder', ...
		'lastUsedRightImageFolder', 'lastUsedStitchedImageFolder', ... 
		'lastUsedOrientationOption', 'lastUsedSizeMismatchOption', ... 
		'lastUsedPauseOption', 'lastUsedOverwriteOption', 'lastUsedSaveFileNameOption');
    return; % SaveParameters


%=====================================================================
% Saves the folders and other user option settings to the mat file.
function handles = RecallParameters(handles)
	clear initialValues; % Need to clear this or else it can have junk in it from a previous running.
	% Load up the initial values from the mat file.
	strIniFile = fullfile(handles.macroFolder, 'Stitch.mat');
	if exist(strIniFile, 'file')
		% Pull out values and stuff them in structure initialValues.
		initialValues = load('Stitch.mat');
		% Assign the image folder from the lastUsedImageFolder field of the
		% structure.
	    handles.leftImageFolder = initialValues.lastUsedLeftImageFolder;
	    handles.rightImageFolder = initialValues.lastUsedRightImageFolder;
	    handles.stitchedImageFolder = initialValues.lastUsedStitchedImageFolder;		
		% Initialize controls with values retrieved from the file.
		set(handles.popOrientation, 'value' , initialValues.lastUsedOrientationOption);
		set(handles.popSizeMismatch, 'value' , initialValues.lastUsedSizeMismatchOption);
		set(handles.chkPauseAfterImage, 'value' , initialValues.lastUsedPauseOption);
		set(handles.chkPromptToOverwrite, 'value' , initialValues.lastUsedOverwriteOption);
		set(handles.chkPromptForFileName, 'value' , initialValues.lastUsedSaveFileNameOption);

		% Default to the current folder if the "last used" one doesn't exist.
		if ~exist(handles.leftImageFolder, 'dir')
			handles.leftImageFolder = cd;
		end
		if ~exist(handles.rightImageFolder, 'dir')
			handles.rightImageFolder = cd;
		end
		if ~exist(handles.stitchedImageFolder, 'dir')
			handles.stitchedImageFolder = cd;
		end
	else
		% If the file is not there, point the image folder to the current
		% directory.
		handles.leftImageFolder = cd;
		handles.rightImageFolder = cd;
		% Make default Stitched folder a sub-folder of the left image folder.
	    handles.stitchedImageFolder = [handles.leftImageFolder '\Stitched'];
	end
    set(handles.txtLeftFolder, 'string' ,handles.leftImageFolder);
    set(handles.txtRightFolder, 'string' ,handles.rightImageFolder);
    set(handles.txtStitchedFolder, 'string' ,handles.stitchedImageFolder);
	
    %uiwait(msgbox(handles.ImageFolder));
    % Load list of images in the left and right image folders.
    handles = LoadImageList(handles, 1, 1, 1);
	% Select none of the items in the listbox.
	set(handles.lstLeftImageList, 'value', []);
	set(handles.lstRightImageList, 'value', []);
	set(handles.lstStitchedImageList, 'value', []);
	% Update the number of images in the Analyze button caption.
	handles = UpdateAnalyzeButtonCaption(handles);
    
	% Make sure the size mismatch options are correct for
	% the selected orientation option.
	SetSizeMismatchPopupEntries(handles);
	
	return; % from RecallParameters


%=====================================================================
% Make sure the size mismatch options are correct for
% the selected orientation option.
function SetSizeMismatchPopupEntries(handles)
	contents = get(handles.popMode,'String'); % returns popup contents as cell array.
	stitchBlendSelection = contents{get(handles.popMode,'Value')};
	if findstr(stitchBlendSelection, 'Blend') > 0
		% User select the "Blend" option.
		set(handles.popOrientation, 'enable', 'off');
		strList = {'Upper left','Centered','Lower Right'};
		set(handles.popSizeMismatch, 'String', strList);
	else
		% User selected the "Stitch" mode.
		set(handles.popOrientation, 'enable', 'on');
		% For the stitch mode, you need to enable the orientation popup
		% and set the "Size mismatch" popup entries depending on what the orientation selection is.
		orientationOption = get(handles.popOrientation, 'value');
		if orientationOption == 1 
			% Side-by-side
			strList = {'Top','Center','Bottom'};
			set(handles.popSizeMismatch, 'String', strList);
		else
			strList = {'Left','Center','Right'};
			set(handles.popSizeMismatch, 'String', strList);
		end
	end
	return; % SetSizeMismatchPopupEntries

%=====================================================================
% --- Load up the listbox with tif files in folder handles.handles.ImageFolder
function handles=LoadImageList(handles, loadLeftListbox, loadRightListbox, loadStitchedListbox)
	if loadLeftListbox == 1
		% Build up left folder.
		ListOfImageNames = {};
		folder = handles.leftImageFolder;
		if ~isempty(folder) && length(folder) > 0 
			if exist(folder,'dir') 
				% If it gets to here, the folder is good.
				ImageFiles = dir([folder '\*.*']);
				for Index = 1:length(ImageFiles)
					baseFileName = ImageFiles(Index).name;
					[folder, name, extension] = fileparts(baseFileName);
					extension = upper(extension);
					switch lower(extension)
					case {'.png', '.bmp', '.jpg', '.tif', '.wav'}
						% Allow only PNG, TIF, JPG, or BMP images
						ListOfImageNames = [ListOfImageNames baseFileName];
					otherwise
					end
				end
			else
				errorMessage = sprintf('Left Folder\n\t%s\ndoes not exist.', folder);
				msgboxw(errorMessage);
			end
		else
			msgboxw('No left folder specified in function LoadImageList.');
		end
		set(handles.lstLeftImageList, 'string', ListOfImageNames);
		set(handles.txtLeftFolder, 'string', handles.leftImageFolder);
		% Need to reset the number of items selected to null otherwise,
		% it will remember the selection from before, even if the list has
		% been cleared out and regenerated.
		set(handles.lstLeftImageList, 'value', []);
	end
	
	if loadRightListbox == 1
		% Do the same for the right folder.
		ListOfImageNames = {};
		folder = handles.rightImageFolder;
		if length(folder) > 0 
			if exist(folder,'dir') 
				% If it gets to here, the folder is good.
				ImageFiles = dir([folder '\*.*']);
				for Index = 1:length(ImageFiles)
					baseFileName = ImageFiles(Index).name;
					[folder, name, extension] = fileparts(baseFileName);
					extension = upper(extension);
					switch lower(extension)
					case {'.png', '.bmp', '.jpg', '.tif', '.wav'}
						% Allow only PNG, TIF, JPG, or BMP images
						ListOfImageNames = [ListOfImageNames baseFileName];
					otherwise
					end
				end
			else
				errorMessage = sprintf('Right Folder\n\t%s\ndoes not exist.', folder);
				msgboxw(errorMessage);
			end
		else
			msgboxw('No right folder specified in function LoadImageList.');
		end
		set(handles.lstRightImageList, 'string', ListOfImageNames);
		set(handles.txtRightFolder, 'string', handles.rightImageFolder);
		% Need to reset the number of items selected to null otherwise,
		% it will remember the selection from before, even if the list has
		% been cleared out and regenerated.
		set(handles.lstRightImageList, 'value', []);
	end
	
	if loadStitchedListbox == 1
		% Do the same for the stitched folder.
		ListOfImageNames = {};
		folder = handles.stitchedImageFolder;
		if length(folder) > 0 
			if exist(folder,'dir') 
				% If it gets to here, the folder is good.
				ImageFiles = dir([folder '\*.*']);
				for Index = 1:length(ImageFiles)
					baseFileName = ImageFiles(Index).name;
					[folder, name, extension] = fileparts(baseFileName);
					extension = upper(extension);
					switch lower(extension)
					case {'.png', '.bmp', '.jpg', '.tif', '.wav'}
						% Allow only PNG, TIF, JPG, or BMP images
						ListOfImageNames = [ListOfImageNames baseFileName];
					otherwise
					end
				end
			else
				errorMessage = sprintf('Stitched Folder\n\t%s\ndoes not exist.', folder);
				msgboxw(errorMessage);
			end
		else
			msgboxw('No stitched folder specified in function LoadImageList.');
		end
		set(handles.lstStitchedImageList, 'string', ListOfImageNames);
		set(handles.txtStitchedFolder, 'string', handles.stitchedImageFolder);
		% Need to reset the number of items selected to null otherwise,
		% it will remember the selection from before, even if the list has
		% been cleared out and regenerated.
		set(handles.lstStitchedImageList, 'value', []);
		% Disable Delete button.
		set(handles.btnDelete, 'Enable', 'off');		
	end
	
    return; % LoadImageList
	

%=====================================================================
% --- Executes on clicking btnStitch button.
% Finds out what they want to do (stitch or blend) and calls the function to do that.
function btnStitch_Callback(hObject, eventdata, handles)
	% Find out if they specified Stitch or Blend.
	selectedIndex = get(handles.popMode, 'Value');
	%mode = contents{selectedIndex}
	if selectedIndex == 1
		% Stitch
		StitchImages(hObject, eventdata, handles);
	else
		% Blend
		BlendImages(hObject, eventdata, handles);
	end

%=====================================================================
% Goes down through the list, displaying then analyzing each highlighted image file.
function StitchImages(hObject, eventdata, handles)
    % hObject    handle to btnStitch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

	% Change mouse pointer (cursor) to an hourglass.  
	% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
	
	global stitchedImage;	% Stitched_Image_Click_Event needs to use this array.

	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.

	% Get a list of those indexes that are selected, so we know which images to process.
    leftSelected = get(handles.lstLeftImageList, 'value');
    numberOfLeftSelectedFiles = uint16(length(leftSelected));
    rightSelected = uint16(get(handles.lstRightImageList, 'value'));
    numberOfRightSelectedFiles = uint16(length(rightSelected));
	
	% There must be the same number of files selected on both sides.
	if numberOfLeftSelectedFiles ~= numberOfRightSelectedFiles
		message = sprintf('You must have the same number of images selected on both sides.\nYou have selected %d images on the left side\nand have selected %d images on the right side', numberOfLeftSelectedFiles, numberOfRightSelectedFiles);
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		msgboxw(message);
		return;
	end

	% Then get list of all of the filenames in the list,
	% regardless of whether they are selected or not.
    ListOfLeftImageNames = get(handles.lstLeftImageList, 'string');
    ListOfRightImageNames = get(handles.lstRightImageList, 'string');

	% Get the orientation that they want over/under or side-by-side.
	orientationOption = get(handles.popOrientation, 'value');
	% side-by-side is 1, over/under is 2.

	if numberOfLeftSelectedFiles >= 2
		% Make progress pie chart visible.
		% Can only make axes disappear if you put it in a panel.
		set(handles.pnlProgress, 'Visible', 'on');
		set(handles.axesPie, 'Visible', 'on');
	end

	imageCount = uint16(0);
	for j=1:numberOfLeftSelectedFiles    % Loop though all selected indexes.
        leftIndex = leftSelected(j);    % Get the next selected index.
        rightIndex = rightSelected(j);    % Get the next selected index.
        % Get the filename for this selected index.
        leftImageName = strcat(cell2mat(ListOfLeftImageNames(leftIndex)));
		leftImageFullFileName = fullfile(handles.leftImageFolder, leftImageName);
		set(handles.txtLeftImageName, 'String', leftImageName);
        rightImageName = strcat(cell2mat(ListOfRightImageNames(rightIndex)));
		rightImageFullFileName = fullfile(handles.rightImageFolder, rightImageName);
		set(handles.txtRightImageName, 'String', rightImageName);
		
		% Get the extensions.  If they're wave files, we need to process a different way, plus they
		% will both need to be wave files to stitch together (you can't stitch a wav file to a jpg
		% file).
		[driveFolderL, BaseNameL, extensionL] = fileparts(leftImageFullFileName) ;
		[driveFolderR, BaseNameR, extensionR] = fileparts(rightImageFullFileName) ;
		if strcmpi(extensionL, '.wav') > 0 || strcmpi(extensionR, '.wav') > 0
			% At least one is a wav file.  If one is, both must be.
			if strcmpi(extensionL, extensionR) <= 0 
				message = sprintf('You cannot append a %s file to a %s file!\n%s\n%s', extensionL, extensionR, leftImageName, rightImageName);
				uiwait(msgbox(message));
			else
				AppendWavFiles(leftImageFullFileName, rightImageFullFileName, handles.stitchedImageFolder);
			end
			continue;	% Skip to end of for loop.
		end
		
		% Read in and display the images (they are not rotated yet).
		leftImage = DisplayImage(handles, 'Left', leftImageFullFileName);
		rightImage = DisplayImage(handles, 'Right', rightImageFullFileName);

		% See how many color planes they have.  We'll convert them both to
		% color later if one is color and one is not color.
		leftImageSize = size(leftImage);
		rightImageSize = size(rightImage);
		if length(leftImageSize) >= 3
			numberOfLeftColors = leftImageSize(3);
		else
			numberOfLeftColors = 1;
		end
		if length(rightImageSize) >= 3
			numberOfRightColors = rightImageSize(3);
		else
			numberOfRightColors = 1;
		end
		
		% Rotate it in the original image thumbnails if the checkbox is checked.
		showRotated = get(handles.chkShowRotated, 'Value');
		
		% Rotate the left image if necessary.
		leftSelectedRotationIndex = get(handles.popLeftRotation,'Value'); % returns selected item from popLeftRotation
		switch(leftSelectedRotationIndex)
			case 2
				% 90 clockwise
				leftImage = imrotate(leftImage, -90);
			case 3
				% 90 counter clockwise
				leftImage = imrotate(leftImage, 90);
			case 4
				% 180
				leftImage = imrotate(leftImage, 180);
		end
		% Display the original image rotated also.
		if showRotated == 1 && leftSelectedRotationIndex ~= 1
			axes(handles.axesLeftImage);
			imshow(leftImage);
		end

		% Rotate the right image if necessary.
		rightSelectedRotationIndex = get(handles.popRightRotation,'Value'); % returns selected item from popRightRotation
		switch(rightSelectedRotationIndex)
			case 2
				% 90 clockwise
				rightImage = imrotate(rightImage, -90);
			case 3
				% 90 counter clockwise
				rightImage = imrotate(rightImage, 90);
			case 4
				% 180
				rightImage = imrotate(rightImage, 180);
		end
		% Display the original image rotated also.
		if showRotated == 1 && rightSelectedRotationIndex ~= 1
			axes(handles.axesRightImage);
			imshow(rightImage);
		end

		% Convert to color if necessary.
		if numberOfRightColors ~= numberOfLeftColors
			% Number of color bands don't match.
			% Let's assume one is 3 (RGB) and the other is 1 (monochrome).
			if numberOfRightColors == 1 
				rightImage = gray2rgb(rightImage);
			elseif numberOfLeftColors == 1
				leftImage = gray2rgb(leftImage);
			end
		end
		
		% Make sure they both have the same number of rows or columns
		% (depends on the orientation.)
		leftImageSize = size(leftImage);
		rightImageSize = size(rightImage);
		if leftImageSize(1) ~= rightImageSize(1) || leftImageSize(2) ~= rightImageSize(2)
			% Dimensions of images don't match.
			expandedCanvas = CreateExpandedCanvass(leftImageSize, rightImageSize, orientationOption);
			expandedCanvasSize = size(expandedCanvas);
			
			% See how user wants to handle it.
			sizeMismatchOption = get(handles.popSizeMismatch, 'Value');
			switch (sizeMismatchOption)
				case 4
					% Skip image but issue no warning.
					continue;
				case 5
					% Skip image with a warning.
					errorMessage = sprintf('The image dimensions do not match.\nLeft image dimensions = %d by %d.\nRight image dimensions = %d by %d', leftImageSize(1), leftImageSize(2), rightImageSize(1), rightImageSize(2));
					msgboxw(errorMessage);
					continue;
			end
			% For cases 1-3, we need to adjust the position of the image
			% to fit into the expanded canvass at the proper position.
			if orientationOption == 1
				% Side-by-side option.
				% Find out which image is taller (has more rows).
				if leftImageSize(1) > rightImageSize(1)
					% Left image is taller.
					switch (sizeMismatchOption)
						case 1
							% Align tops of images.
							topLine = 1;
							bottomLine = topLine + rightImageSize(1) - 1;
						case 2
							% Align centers of images.
							topLine = floor((leftImageSize(1) - rightImageSize(1)) / 2);
							bottomLine = topLine + rightImageSize(1) - 1;
						case 3
							% Align bottoms of images.
							topLine = leftImageSize(1) - rightImageSize(1) + 1;
							bottomLine = topLine + rightImageSize(1) - 1;
						otherwise
							% Shouldn't happen.
							continue;
					end
					% Insert rightImage into the expanded canvass.
					% They should both be either color or monochrome at this point.
					if length(rightImageSize) <= 2
						% Monochrome image.
						expandedCanvas(topLine:bottomLine,1:rightImageSize(2)) = rightImage;
					else
						% Color image.
						expandedCanvas(topLine:bottomLine,1:rightImageSize(2), :) = rightImage;
					end
					% Replace the expanded canvass with the new, larger rightImage.
					rightImage = expandedCanvas;
				else
					% Right image is taller.
					switch (sizeMismatchOption)
						case 1
							% Align tops of images.
							topLine = 1;
							bottomLine = topLine + leftImageSize(1) - 1;
						case 2
							% Align centers of images.
							topLine = floor((rightImageSize(1) - leftImageSize(1)) / 2);
							if topLine <= 0; topLine = 1; end;
							bottomLine = topLine + leftImageSize(1) - 1;
						case 3
							% Align bottoms of images.
							topLine = rightImageSize(1) - leftImageSize(1) + 1;
							if topLine <= 0; topLine = 1; end;
							bottomLine = topLine + leftImageSize(1) - 1;
						otherwise
							% Shouldn't happen.
							continue;
					end
					% Insert leftImage into the expanded canvass.
					if length(rightImageSize) <= 2
						% Monochrome image.
						expandedCanvas(topLine:bottomLine,1:leftImageSize(2)) = leftImage;
					else
						% Color image.
						expandedCanvas(topLine:bottomLine,1:leftImageSize(2), :) = leftImage;
					end
					% Replace the expanded canvass with the new, larger leftImage.
					leftImage = expandedCanvas;
				end
				
			else
				% Over/under option.
				topLine = 1;
				bottomLine = expandedCanvasSize(1);
				% Find out which image is wider (has more columns).
				if leftImageSize(2) > rightImageSize(2)
					% Left (top) image is wider.
					switch (sizeMismatchOption)
						case 1
							% Align lefts of images.
							leftColumn = 1;
							rightColumn = leftColumn + rightImageSize(2) - 1;
						case 2
							% Align centers of images.
							leftColumn = floor((leftImageSize(2) - rightImageSize(2)) / 2);
							if leftColumn <= 0; leftColumn = 1; end;
							rightColumn = leftColumn + rightImageSize(2) - 1;
						case 3
							% Align rights of images.
							leftColumn = leftImageSize(2) - rightImageSize(2) + 1;
							if leftColumn <= 0; leftColumn = 1; end;
							rightColumn = leftColumn + rightImageSize(2) - 1;
						otherwise
							% Shouldn't happen.
							continue;
					end
					% Insert rightImage into the expanded canvass.
					% They should both be either color or monochrome at this point.
					if length(rightImageSize) <= 2
						% Monochrome image.
						expandedCanvas(topLine:bottomLine, leftColumn:rightColumn) = rightImage;
					else
						% Color image.
						expandedCanvas(topLine:bottomLine, leftColumn:rightColumn, :) = rightImage;
					end
					% Replace the expanded canvass with the new, larger rightImage.
					rightImage = expandedCanvas;
				else
					% Right image is wider.
					switch (sizeMismatchOption)
						case 1
							% Align lefts of images.
							leftColumn = 1;
							rightColumn = leftColumn + leftImageSize(2) - 1;
						case 2
							% Align centers of images.
							leftColumn = floor((rightImageSize(2) - leftImageSize(2)) / 2);
							if leftColumn <= 0; leftColumn = 1; end;
							rightColumn = leftColumn + leftImageSize(2) - 1;
						case 3
							% Align rights of images.
							leftColumn = rightImageSize(2) - leftImageSize(2) + 1;
							if leftColumn <= 0; leftColumn = 1; end;
							rightColumn = leftColumn + leftImageSize(2) - 1;
						otherwise
							% Shouldn't happen.
							continue;
					end
					% Insert leftImage into the expanded canvass.
					if length(rightImageSize) <= 2
						% Monochrome image.
						expandedCanvas(topLine:bottomLine, leftColumn:rightColumn) = leftImage;
					else
						% Color image.
						expandedCanvas(topLine:bottomLine, leftColumn:rightColumn, :) = leftImage;
					end
					% Replace the expanded canvass with the new, larger leftImage.
					leftImage = expandedCanvas;
				end
			end
		end

		% Do the actual stitching together
		if orientationOption == 1
			% Side-by-side option.
			stitchedImage = [leftImage rightImage];
		else
			% Over/under option.
			stitchedImage = [leftImage; rightImage];
		end

		% Put stitched image into axesStitchedImage.
		axes(handles.axesStitchedImage);
		handleToImage = imshow(stitchedImage);
		% Make it so that when the click on the image, it will come up in imtool.
		set(handleToImage, 'ButtonDownFcn', @Stitched_Image_Click_Event);
		
		% Save the image in the stitched folder.
		promptForFileName = get(handles.chkPromptForFileName, 'Value');
		userCanceled = 0;
		if promptForFileName
			% Ask user for filename.
			[baseFileName, folder] = uiputfile([handles.stitchedImageFolder '\*.*'],'Save file name');
			if baseFileName ~= 0
				% User didn't cancel out.
				handles.stitchedImageFolder = folder;
				% Use same name as the left image.
				strOutputFullFileName = fullfile(folder, baseFileName);			
			else
				userCanceled = 1;	% Flag to skip writing the file.
			end
		else
			% Use a name like "LeftImageName_plus_RightImageName.ext".
			stitchedImageName = [BaseNameL '_plus_' BaseNameR extensionL];
			strOutputFullFileName = fullfile(handles.stitchedImageFolder, stitchedImageName);			
		end
		status = 1;	% Assume directory exists.
		if ~exist(handles.stitchedImageFolder, 'dir')  && userCanceled == 0
			% Folder doesn't exist yet.  Need to make it.
			[status, message, messageID] = mkdir(handles.stitchedImageFolder);
		end
		if status == 1 && userCanceled == 0
			% Write image file if the folder already exists or 
			% we were successful in creating the new folder.
			% First see if we need to prompt to overwrite.
			if exist(strOutputFullFileName, 'file') && promptForFileName == 0
				% File exists already, but we weren't yet asked about overwriting.
				% If Prompt for file name was checked, they've already been
				% asked about overwriting by the time they get here.
				promptToOverwrite = get(handles.chkPromptToOverwrite, 'value');
				if promptToOverwrite
					% Exists and they want to be prompted.
					promptString = sprintf('File already exists:\n %s\nOverwrite it?', strOutputFullFileName);
					button = questdlg(promptString,'Overwrite?','Yes','No','Yes');
					if strcmpi(button, 'Yes')
						imwrite(stitchedImage, strOutputFullFileName);
					end
				else
					% Exists but they don't want to be prompted.
					imwrite(stitchedImage, strOutputFullFileName);
				end
			else
				% File doesn't exist yet.  Create it.
				imwrite(stitchedImage, strOutputFullFileName);
			end
			set(handles.txtStitchedImageName, 'String', strOutputFullFileName);		
		end

		% Show a pie chart as a progress indicator.
		if numberOfLeftSelectedFiles >= 2
			fractionDone = double(j) / double(numberOfLeftSelectedFiles);
			axes(handles.axesPie);
			pie([fractionDone, (1-fractionDone)]);
		end
		imageCount = imageCount + 1;
		% Prompt to allow user to inspect the image.
		if j < numberOfLeftSelectedFiles && get(handles.chkPauseAfterImage, 'value') 
			msgboxw({'Check out results, then', 'click ok to process the next image.'});    
		end
	end
	
	% Update the stitched image listbox.
	handles = LoadImageList(handles, 0, 0, 1);

	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
	
	% Let them know it's done if there's more than 1 that was processed.
	if imageCount > 1
		doneMessage = sprintf('All done stitching together %d pairs of images.', imageCount);
		msgboxw(doneMessage);
	end

	% Make progress pie chart invisible.
	% DOESN'T SEEM TO WORK.
	set(handles.axesPie, 'Visible', 'off');
	% Can only make axes disappear if you put it in a panel.
	set(handles.pnlProgress, 'Visible', 'off');

	% Clear out & clean up some.
	clear('expandedCanvas', 'leftImage', 'rightImage');

	guidata(hObject, handles);
	return  % Stitch

%=====================================================================
% Append two wave files and create a file for the output.
function AppendWavFiles(wavFileName1, wavFileName2, outputFolder)
	% They're both wave files.
	[driveFolderL, BaseNameL, extensionL] = fileparts(wavFileName1) ;
	[driveFolderR, BaseNameR, extensionR] = fileparts(wavFileName2) ;
	if strcmpi(extensionL, '.wav') > 0 || strcmpi(extensionR, '.wav') > 0
		% At least one is a wav file.  If one is, both must be.
		if strcmpi(extensionL, extensionR) <= 0 
			message = sprintf('You cannot append a %s file to a %s file!\n%s\n%s', extensionL, extensionR, wavFileName1, wavFileName2);
			uiwait(msgbox(message));
		else
			% They're both wav files.
			[wavData1, Fs1, nBits1] = wavread(wavFileName1);
			dimensions1 = size(wavData1);
			[wavData2, Fs2, nBits2]  = wavread(wavFileName2);
			dimensions2 = size(wavData2);
			% Make sure if one is stereo, they're both stereo.  You can't append otherwise
			% because the arrays are different dimensions.
			if dimensions1(2) == 2 && dimensions2(2) == 1 
				% 1 is stereo and 2 is mono.
				% Convert wavData2 to stereo.
				stereo2 = [wavData2 wavData2];
				both = [wavData1; stereo2];
			elseif dimensions1(2) == 1 && dimensions2(2) == 2
				% 1 is stereo and 2 is mono.
				% Convert wavData1 to stereo.
				stereo1 = [wavData1 wavData1];
				both = [stereo1; wavData2];
			else
				% They're both the same - either both stereo or both mono.
				both = [wavData1; wavData2];
			end
			% Use same name as the left image.
			stitchedName = sprintf('%s_plus_%s.wav',BaseNameL, BaseNameR);
			strOutputFullFileName = fullfile(outputFolder, stitchedName);			
			wavwrite(both, Fs1, strOutputFullFileName);
		end
	end
	return  % AppendWavFiles
	
%=====================================================================
% Goes down through the list, displaying then analyzing each highlighted image file.
function BlendImages(hObject, eventdata, handles)
    % hObject    handle to btnStitch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

	% Change mouse pointer (cursor) to an hourglass.  
	% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
	
	global stitchedImage;	% Stitched_Image_Click_Event needs to use this array.

	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.

	% Get a list of those indexes that are selected, so we know which images to process.
    leftSelected = get(handles.lstLeftImageList, 'value');
    numberOfLeftSelectedFiles = uint16(length(leftSelected));
    rightSelected = uint16(get(handles.lstRightImageList, 'value'));
    numberOfRightSelectedFiles = uint16(length(rightSelected));
	
	% There must be the same number of files selected on both sides.
	if numberOfLeftSelectedFiles ~= numberOfRightSelectedFiles
		message = sprintf('You must have the same number of images selected on both sides.\nYou have selected %d images on the left side\nand have selected %d images on the right side', numberOfLeftSelectedFiles, numberOfRightSelectedFiles);
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		msgboxw(message);
		return;
	end

	% Then get list of all of the filenames in the list,
	% regardless of whether they are selected or not.
    ListOfLeftImageNames = get(handles.lstLeftImageList, 'string');
    ListOfRightImageNames = get(handles.lstRightImageList, 'string');

	% Get the orientation that they want over/under or side-by-side.
	orientationOption = get(handles.popOrientation, 'value');
	% side-by-side is 1, over/under is 2.

	if numberOfLeftSelectedFiles >= 2
		% Make progress pie chart visible.
		% Can only make axes disappear if you put it in a panel.
		set(handles.pnlProgress, 'Visible', 'on');
		set(handles.axesPie, 'Visible', 'on');
	end

	imageCount = uint16(0);
	for j = 1 : numberOfLeftSelectedFiles    % Loop though all selected indexes.
        leftIndex = leftSelected(j);    % Get the next selected index.
        rightIndex = rightSelected(j);    % Get the next selected index.
        % Get the filename for this selected index.
        leftImageName = strcat(cell2mat(ListOfLeftImageNames(leftIndex)));
		leftImageFullFileName = fullfile(handles.leftImageFolder, leftImageName);
		set(handles.txtLeftImageName, 'String', leftImageName);
        rightImageName = strcat(cell2mat(ListOfRightImageNames(rightIndex)));
		rightImageFullFileName = fullfile(handles.rightImageFolder, rightImageName);
		set(handles.txtRightImageName, 'String', rightImageName);
		
		% I don't have the ability to blend sound files at this time, 
		% so bail out if they selected any wav files.
		waveFileOnLeft = findstr(leftImageName, '.wav');
		waveFileOnRight = findstr(rightImageName, '.wav');
		% The above will be empty [] unless there's .wav in the filename.
		if ~isempty(waveFileOnLeft) || ~isempty(waveFileOnRight)
			notice = sprintf('This program does not have the ability to blend wave files.\nYou can only stitch them.');
			msgboxw(notice);
			set(gcf,'Pointer','arrow');
			drawnow;	% Cursor won't change right away unless you do this.
			return;
		end

		% Read in and display the images (they are not rotated yet).
		leftImage = DisplayImage(handles, 'Left', leftImageFullFileName);
		rightImage = DisplayImage(handles, 'Right', rightImageFullFileName);
		
		% See how many color planes they have.  We'll convert them both to
		% color later if one is color and one is not color.
		leftImageSize = size(leftImage);
		rightImageSize = size(rightImage);
		if length(leftImageSize) >= 3
			numberOfLeftColors = leftImageSize(3);
		else
			numberOfLeftColors = 1;
		end
		if length(rightImageSize) >= 3
			numberOfRightColors = rightImageSize(3);
		else
			numberOfRightColors = 1;
		end
		
		% Rotate it in the original image thumbnails if the checkbox is checked.
		showRotated = get(handles.chkShowRotated, 'Value');
		
		% Rotate the left image if necessary.
		leftSelectedRotationIndex = get(handles.popLeftRotation,'Value'); % returns selected item from popLeftRotation
		switch(leftSelectedRotationIndex)
			case 2
				% 90 clockwise
				leftImage = imrotate(leftImage, -90);
			case 3
				% 90 counter clockwise
				leftImage = imrotate(leftImage, 90);
			case 4
				% 180
				leftImage = imrotate(leftImage, 180);
		end
		% Display the original image rotated also.
		if showRotated == 1 && leftSelectedRotationIndex ~= 1
			axes(handles.axesLeftImage);
			imshow(leftImage);
		end

		% Rotate the right image if necessary.
		rightSelectedRotationIndex = get(handles.popRightRotation,'Value'); % returns selected item from popRightRotation
		switch(rightSelectedRotationIndex)
			case 2
				% 90 clockwise
				rightImage = imrotate(rightImage, -90);
			case 3
				% 90 counter clockwise
				rightImage = imrotate(rightImage, 90);
			case 4
				% 180
				rightImage = imrotate(rightImage, 180);
		end
		% Display the original image rotated also.
		if showRotated == 1 && rightSelectedRotationIndex ~= 1
			axes(handles.axesRightImage);
			imshow(rightImage);
		end

		% Convert to color if necessary.
		if numberOfRightColors ~= numberOfLeftColors
			% Number of color bands don't match.
			% Let's assume one is 3 (RGB) and the other is 1 (monochrome).
			if numberOfRightColors == 1 
				rightImage = gray2rgb(rightImage);
			elseif numberOfLeftColors == 1
				leftImage = gray2rgb(leftImage);
			end
		end
		
		% Create an image that can fit the maximum dimensions of both images.
		leftImageSize = size(leftImage);
		rightImageSize = size(rightImage);
		if length(leftImageSize) <= 2 && length(rightImageSize) <= 2
			% Both images are monochrome.
			expandedCanvas = zeros(max([leftImageSize(1) rightImageSize(1)]), max([leftImageSize(2) rightImageSize(2)]));
		else
			% At least one image is color, so the output must also be made color.
			expandedCanvas = zeros(max([leftImageSize(1) rightImageSize(1)]), max([leftImageSize(2) rightImageSize(2)]), 3);
		end
		expandedCanvasSize = size(expandedCanvas);

		% See how user wants to handle it.
		sizeMismatchOption = get(handles.popSizeMismatch, 'Value');
		switch (sizeMismatchOption)
			case 4
				% Skip image but issue no warning.
				continue;
			case 5
				% Skip image with a warning.
				errorMessage = sprintf('The image dimensions do not match.\nLeft image dimensions = %d by %d.\nRight image dimensions = %d by %d', leftImageSize(1), leftImageSize(2), rightImageSize(1), rightImageSize(2));
				msgboxw(errorMessage);
				continue;
		end
		% Specify left column and top line for each image.
		switch (sizeMismatchOption)
			case 1
				% Align top lefts of images.
				topLine1 = 1;
				topLine2 = 1;
				leftColumn1 = 1;
				leftColumn2 = 1;
			case 2
				topLine1 = floor((expandedCanvasSize(1) - leftImageSize(1)) / 2);
				leftColumn1 = floor((expandedCanvasSize(2) - leftImageSize(2)) / 2);
				topLine2 = floor((expandedCanvasSize(1) - rightImageSize(1)) / 2);
				leftColumn2 = floor((expandedCanvasSize(2) - rightImageSize(2)) / 2);
			case 3
				topLine1 = expandedCanvasSize(1) - leftImageSize(1) + 1;
				leftColumn1 = expandedCanvasSize(2) - leftImageSize(2) + 1;
				topLine2 = expandedCanvasSize(1) - rightImageSize(1) + 1;
				leftColumn2 = expandedCanvasSize(2) - rightImageSize(2) + 1;
			otherwise
				% Shouldn't happen.
				continue;
		end
		% Make sure none are less than 1 - if so, make them 1.
		[topLine1 topLine2 leftColumn1 leftColumn2] = ClipToOne(topLine1, topLine2, leftColumn1, leftColumn2);
		rightColumn1 = leftColumn1 + leftImageSize(2) - 1;
		rightColumn2 = leftColumn2 + rightImageSize(2) - 1;
		bottomLine1 = topLine1 + leftImageSize(1) - 1;
		bottomLine2 = topLine2 + rightImageSize(1) - 1;
		% Make sure none are less than 1 or greater than the canvass size.
		[topLine1 topLine2 bottomLine1 bottomLine2 leftColumn1 leftColumn2 rightColumn1 rightColumn2] = ...
			CheckBounds(expandedCanvasSize, topLine1, topLine2, bottomLine1, bottomLine2, leftColumn1, leftColumn2, rightColumn1, rightColumn2);
		
		% Insert leftImage into the expanded canvass.
		if length(leftImageSize) <= 2 && length(rightImageSize) <= 2
			% Monochrome image.
			monoImage = single(expandedCanvas);
			% Assign in the first image.
			monoImage(topLine1:bottomLine1, leftColumn1:rightColumn1) = leftImage(:,:);
			% Add in the second image.
			monoImageR = single(rightImage(:,:));
			monoImage(topLine2:bottomLine2, leftColumn2:rightColumn2) = ...
				monoImage(topLine2:bottomLine2, leftColumn2:rightColumn2) + monoImageR;
			% Divide by 2 to get the average.
			expandedCanvas = monoImage / 2.0;
			clear('monoImage', 'monoImageR');
		else
			% Color image.
			% The following line doesn't work because of the following error:
			% ??? Assignment has more non-singleton rhs dimensions than non-singleton subscripts.
			% So it's commented out and we'll do it band by band instead.
% 					expandedCanvas(topLine:bottomLine,1:leftImageSize(2), :) = leftImage; % Doesn't work.
			redBand = single(expandedCanvas(:,:,1));
			greenBand = single(expandedCanvas(:,:,2));
			blueBand = single(expandedCanvas(:,:,3));
			% Assign in the first image.
			redBand(topLine1:bottomLine1, leftColumn1:rightColumn1) = leftImage(:,:,1);
			greenBand(topLine1:bottomLine1, leftColumn1:rightColumn1) = leftImage(:,:,2);
			blueBand(topLine1:bottomLine1, leftColumn1:rightColumn1) = leftImage(:,:,3);
			% Add in the second image.
			redBandR = single(rightImage(:,:,1));
			greenBandR = single(rightImage(:,:,2));
			blueBandR = single(rightImage(:,:,3));
			redBand(topLine2:bottomLine2, leftColumn2:rightColumn2) = ...
				redBand(topLine2:bottomLine2, leftColumn2:rightColumn2) + redBandR;
			greenBand(topLine2:bottomLine2, leftColumn2:rightColumn2) = ...
				greenBand(topLine2:bottomLine2, leftColumn2:rightColumn2) + greenBandR;
			blueBand(topLine2:bottomLine2, leftColumn2:rightColumn2) = ...
				blueBand(topLine2:bottomLine2, leftColumn2:rightColumn2) + blueBandR;
			% Divide by 2 to get the average.
			redBand = redBand / 2.0;
			greenBand = greenBand / 2.0;
			blueBand = blueBand / 2.0;
			expandedCanvas = cat(3, redBand, greenBand, blueBand);
			clear('redBand', 'greenBand', 'blueBand');
		end
		% Now the expanded canvass should have the average of both images.
		stitchedImage = uint8(expandedCanvas);

		% Put blended image into axesStitchedImage.
		axes(handles.axesStitchedImage);
		handleToImage = imshow(stitchedImage);
		% Make it so that when the click on the image, it will come up in imtool.
		set(handleToImage, 'ButtonDownFcn', @Stitched_Image_Click_Event);
		
		% Save the image in the stitched folder.
		promptForFileName = get(handles.chkPromptForFileName, 'Value');
		userCanceled = 0;
		if promptForFileName
			% Ask user for filename.
			[baseFileName, folder] = uiputfile([handles.stitchedImageFolder '\*.*'],'Save file name');
			if baseFileName ~= 0
				% User didn't cancel out.
				handles.stitchedImageFolder = folder;
				% Use same name as the left image.
				strOutputFullFileName = fullfile(folder, baseFileName);			
			else
				userCanceled = 1;	% Flag to skip writing the file.
			end
		else
			% Create a new filename from the left and right image names.
			[folder, leftBaseFileName, extension] = fileparts(leftImageName);
			[folder, rightBaseFileName, extension] = fileparts(rightImageName);
			strOutputFullFileName = sprintf('Mean_of_%s_and_%s.tif', leftBaseFileName, rightBaseFileName);
			strOutputFullFileName = fullfile(handles.stitchedImageFolder, strOutputFullFileName);			
		end
		status = 1;	% Assume directory exists.
		if ~exist(handles.stitchedImageFolder, 'dir')  && userCanceled == 0
			% Folder doesn't exist yet.  Need to make it.
			[status, message, messageID] = mkdir(handles.stitchedImageFolder);
		end
		if status == 1 && userCanceled == 0
			% Write image file if the folder already exists or 
			% we were successful in creating the new folder.
			% First see if we need to prompt to overwrite.
			if exist(strOutputFullFileName, 'file') && promptForFileName == 0
				% File exists already, but we weren't yet asked about overwriting.
				% If Prompt for file name was checked, they've already been
				% asked about overwriting by the time they get here.
				promptToOverwrite = get(handles.chkPromptToOverwrite, 'value');
				if promptToOverwrite
					% Exists and they want to be prompted.
					promptString = sprintf('File already exists:\n %s\nOverwrite it?', strOutputFullFileName);
					button = questdlg(promptString,'Overwrite?','Yes','No','Yes');
					if strcmpi(button, 'Yes')
						imwrite(stitchedImage, strOutputFullFileName);
					end
				else
					% Exists but they don't want to be prompted.
					imwrite(stitchedImage, strOutputFullFileName);
				end
			else
				% File doesn't exist yet.  Create it.
				imwrite(stitchedImage, strOutputFullFileName);
			end
			set(handles.txtStitchedImageName, 'String', strOutputFullFileName);		
		end

		% Show a pie chart as a progress indicator.
		if numberOfLeftSelectedFiles >= 2
			fractionDone = double(j) / double(numberOfLeftSelectedFiles);
			axes(handles.axesPie);
			pie([fractionDone, (1-fractionDone)]);
		end
		imageCount = imageCount + 1;
		% Prompt to allow user to inspect the image.
		if j < numberOfLeftSelectedFiles && get(handles.chkPauseAfterImage, 'value') 
			msgboxw({'Check out results, then', 'click ok to process the next image.'});    
		end
	end
	
	% Update the stitched image listbox.
	handles = LoadImageList(handles, 0, 0, 1);

	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
	
	% Let them know it's done if there's more than 1 that was processed.
	if imageCount > 1
		doneMessage = sprintf('All done blending together %d pairs of images.', imageCount);
		msgboxw(doneMessage);
	end

	% Make progress pie chart invisible.
	% DOESN'T SEEM TO WORK.
	set(handles.axesPie, 'Visible', 'off');
	% Can only make axes disappear if you put it in a panel.
	set(handles.pnlProgress, 'Visible', 'off');

	% Clear out & clean up some.
	clear('expandedCanvas', 'leftImage', 'rightImage');
	guidata(hObject, handles);
	return;  % BlendImages

%=====================================================================
% Creates an expanded canvass that is large enough to hold both images.
% Correctly handles both side-by-side option and over/under option.
% Canvass is simply initialized to all zeros and will later get elements
% replaced in the calling routine.
function expandedCanvas = CreateExpandedCanvass(leftImageSize, rightImageSize, orientationOption)
	if orientationOption == 1
		% Side-by-side option.
		% Find out which image has more rows.
		if leftImageSize(1) > rightImageSize(1)
			% Left image is taller.
			% Canvass needs to get taller to hold the right image,
			% but retains the same width as the right image already has.
			expandedCanvassHeight = leftImageSize(1);
			expandedCanvassWidth = rightImageSize(2);
		else
			% Right image is taller.
			% Canvass needs to get taller to hold the left image,
			% but retains the same width as the left image already has.
			expandedCanvassHeight = rightImageSize(1);
			expandedCanvassWidth = leftImageSize(2);
		end
	else
		% Over/under option.
		% Find out which image has more columns.
		if leftImageSize(2) > rightImageSize(2)
			% Left image is wider.
			% Canvass needs to get wider to hold the right image,
			% but retains the same height as the right image already has.
			expandedCanvassHeight = rightImageSize(1);
			expandedCanvassWidth = leftImageSize(2);
		else
			% Right image is wider.
			% Canvass needs to get wider to hold the left image,
			% but retains the same height as the left image already has.
			expandedCanvassHeight = leftImageSize(1);
			expandedCanvassWidth = rightImageSize(2);
		end
	end
	% Make a new image with all zeros (for the time being) that
	% can fit both images.  If at least one of the images is color, make
	% the canvass color, otherwise make it monochrome.
	if length(leftImageSize) > 2 || length(rightImageSize) > 2
		% At least one image is a color image - make a 3D array.
		expandedCanvas = zeros(expandedCanvassHeight, expandedCanvassWidth, 3);
	else
		% They're both monochrome images - make a 2D array.
		expandedCanvas = zeros(expandedCanvassHeight, expandedCanvassWidth);
	end
	return; % CreateExpandedCanvass
	
%=====================================================================
% Converts a 2D, 1 plane image grayImage(:,:) into a
% 2D, 3 plane image colorImage(:,:,3)
function colorImage = gray2rgb(grayImage)
	colorImage = repmat(grayImage,[1,1,3]);
	return;
		
%=====================================================================
% --- Executes on button press in btnSelectAllOrNone.
function btnSelectAllOrNone_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectAllOrNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Find out button caption and take appropriate action.
    ButtonCaption = get(handles.btnSelectAllOrNone, 'string');
	if strcmp(ButtonCaption, 'Step 2:  Select All') == 1
        % Select all items in the listbox.
        % Need to find out how many items are in the listbox (both selected and
        % unselected).  It's quirky and inefficient but it's the only way I
        % know how to do it. 
        % First get the whole damn listbox text into a cell array.
        caListboxString = get(handles.lstLeftImageList, 'string');
        NumberOfItems = length(caListboxString);    % Get length of that cell array.
        AllIndices=1:NumberOfItems; % Make a vector of all indices.
        % Select all indices.
        set(handles.lstLeftImageList, 'value', AllIndices);
        % First get the whole damn listbox text into a cell array.
        caListboxString = get(handles.lstRightImageList, 'string');
        NumberOfItems = length(caListboxString);    % Get length of that cell array.
        AllIndices=1:NumberOfItems; % Make a vector of all indices.
        % Select all indices.
        set(handles.lstRightImageList, 'value', AllIndices);
        % Finally, change caption to say "Select None"
        set(handles.btnSelectAllOrNone, 'string', 'Step 2:  Select None');
        % It scrolls to the bottom of the list.  Use the following line
        % if you want the first item at the top of the list.
        set(handles.lstLeftImageList, 'ListboxTop', 1);
        set(handles.lstRightImageList, 'ListboxTop', 1);
    else
        % Select none of the items in the listbox.
        set(handles.lstLeftImageList, 'value', []);
        set(handles.lstRightImageList, 'value', []);
        % Change caption to say Select All
        set(handles.btnSelectAllOrNone, 'string', 'Step 2:  Select All');
	end
	% Update the number of images in the Analyze button caption.
	handles = UpdateAnalyzeButtonCaption(handles);
	guidata(hObject, handles);


%=====================================================================
% --- Executes on button press in chkPauseAfterImage.
function chkPauseAfterImage_Callback(hObject, eventdata, handles)
% hObject    handle to chkPauseAfterImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPauseAfterImage


%=====================================================================
% --- Executes during object creation, after setting all properties.
% EVEN THOUGH THIS FUNCTION IS EMPTY, DON'T DELETE IT OR ERRORS WILL OCCUR
function axesStitchedImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axesStitchedImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesStitchedImage


%=====================================================================
% --- Executes during object creation, after setting all properties.
% EVEN THOUGH THIS FUNCTION IS EMPTY, DON'T DELETE IT OR ERRORS WILL OCCUR
function axesLeftImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axesLeftImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesLeftImage




%=====================================================================
% --- Executes during object creation, after setting all properties.
% EVEN THOUGH THIS FUNCTION IS EMPTY, DON'T DELETE IT OR ERRORS WILL OCCUR
function figMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%=====================================================================
% --- Executes on button press in btnExit.
function btnExit_Callback(hObject, eventdata, handles)
% hObject    handle to btnExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	% Save all the current settings.
	handles = SaveParameters(handles);
	delete(handles.figMain);


% --- Executes on button press in chkAutoAnalyze.
function chkAutoAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to chkAutoAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkAutoAnalyze


% --- Executes on selection change in popLeftRotation.
function popLeftRotation_Callback(hObject, eventdata, handles)
% hObject    handle to popLeftRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popLeftRotation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLeftRotation


% --- Executes during object creation, after setting all properties.
function popLeftRotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLeftRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popRightRotation.
function popRightRotation_Callback(hObject, eventdata, handles)
% hObject    handle to popRightRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popRightRotation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popRightRotation


% --- Executes during object creation, after setting all properties.
function popRightRotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popRightRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function lstRightImageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstRightImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radRightImageFolder.
function radRightImageFolder_Callback(hObject, eventdata, handles)
% hObject    handle to radRightImageFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radRightImageFolder
    radioButtonSetting = get(hObject,'Value');
    % Set the other radio button's value the opposite of this one.
    if radioButtonSetting == 1
		set(handles.btnSelectFolder, 'String', 'Select Right Image Folder...');
    end
    guidata(hObject, handles);


% --- Executes on button press in radLeftImageFolder.
function radLeftImageFolder_Callback(hObject, eventdata, handles)
% hObject    handle to radLeftImageFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radLeftImageFolder
    radioButtonSetting = get(hObject,'Value');
    % Set the other radio button's value the opposite of this one.
    if radioButtonSetting == 1
%         set(handles.radRightImageFolder, 'value', 0);
		set(handles.btnSelectFolder, 'String', 'Select Left Image Folder...');
    end
    guidata(hObject, handles);


% --- Executes on button press in radStitchedImageFolder.
function radStitchedImageFolder_Callback(hObject, eventdata, handles)
% hObject    handle to radStitchedImageFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radStitchedImageFolder
    radioButtonSetting = get(hObject,'Value');
    % Set the other radio button's value the opposite of this one.
    if radioButtonSetting == 1
		set(handles.btnSelectFolder, 'String', 'Select Stitched Image Folder...');
    end
    guidata(hObject, handles);



% --- Executes on selection change in lstStitchedImageList.
function lstStitchedImageList_Callback(hObject, eventdata, handles)
% hObject    handle to lstStitchedImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns lstStitchedImageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstStitchedImageList
	axes(handles.axesStitchedImage);

	% Change mouse pointer (cursor) to an hourglass.  
	% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.

	% Get list of selected items.
    Selected = get(handles.lstStitchedImageList, 'value');

	if length(Selected) == 0
		% None selected.
		% Disable Delete button.
		set(handles.btnDelete, 'Enable', 'off');
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		guidata(hObject, handles);
		return; % lstStitchedImageList_Callback
	else
		% Enable Delete button.
		set(handles.btnDelete, 'Enable', 'on');		
	end

	% If more than one is selected, bail out.
	if length(Selected) > 1 
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow')
		drawnow;	% Cursor won't change right away unless you do this.
        return;
	end

	% If only one is selected, display it.
    ListOfImageNames = get(handles.lstStitchedImageList, 'string');
	baseImageFileName = strcat(cell2mat(ListOfImageNames(Selected)));
    fullImageFileName = [handles.stitchedImageFolder '\' baseImageFileName];	% Prepend folder.
	set(handles.txtStitchedImageName, 'String', baseImageFileName);
	
	[folder, baseFileName, extension] = fileparts(fullImageFileName);
	switch lower(extension)
	case '.wav'
		% We're playing a wave file, not showing an image.
		try
			% Hide the prior image.
% 			set(handles.axesStitchedImage, 'Visible', 'off');   % Doesn't work!
% 			childHandles = get(handles.axesStitchedImage, 'children');
% 			delete(childHandles);
			[y, Fs, nbits] = wavread(fullImageFileName);
			plot(y);
			wavplay(y,Fs, 'async');
		catch
			strError = lasterror;
			strErrorMessage = sprintf('MATLAB does not support this type of wave file.\n\n%s', strError.message);
			msgboxw(strErrorMessage);
		end
	case {'.mov', '.wmv', '.asf', '.avi'}
		msgboxw('Video files are not supported by Stitch.');
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		return;
	otherwise
		% Display the image.
		stitchedImage = DisplayImage(handles, 'Stitched', fullImageFileName);
	end
	     
	% Change mouse pointer (cursor) to an arrow.
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
    guidata(hObject, handles);
    return; % lstStitchedImageList_Callback



% --- Executes during object creation, after setting all properties.
function lstStitchedImageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstStitchedImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chkShowRotated.
function chkShowRotated_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowRotated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowRotated




% --- Executes on selection change in popSizeMismatch.
function popSizeMismatch_Callback(hObject, eventdata, handles)
% hObject    handle to popSizeMismatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popSizeMismatch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popSizeMismatch


% --- Executes during object creation, after setting all properties.
function popSizeMismatch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSizeMismatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popOrientation.
function popOrientation_Callback(hObject, eventdata, handles)
% hObject    handle to popOrientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popOrientation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popOrientation
	% Get the orientation that they want: over/under or side-by-side.
	orientationOption = get(handles.popOrientation, 'value');
	if orientationOption == 1 
		% Side-by-side
		strList = {'Top','Center','Bottom'};
		set(handles.popSizeMismatch, 'String', strList);
	else
		strList = {'Left','Center','Right'};
		set(handles.popSizeMismatch, 'String', strList);
	end
	


% --- Executes during object creation, after setting all properties.
function popOrientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popOrientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSwapLeftAndRight.
function btnSwapLeftAndRight_Callback(hObject, eventdata, handles)
% hObject    handle to btnSwapLeftAndRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	% Select none of the items in the listbox.
	% Otherwise selected index stays the same and if one list is shorter
	% than the other list, you'll get an error.
	set(handles.lstLeftImageList, 'value', []);
	set(handles.lstRightImageList, 'value', []);
	set(handles.lstStitchedImageList, 'value', []);

	% Swap folder names
	leftFolder = handles.leftImageFolder;
	handles.leftImageFolder = handles.rightImageFolder;
	handles.rightImageFolder = leftFolder;
	
	% Update the left and right image listboxes.
	handles = LoadImageList(handles, 1, 1, 0);

	% Save the image folders in our initialization mat file.
	handles = SaveParameters(handles);
    guidata(hObject, handles);
	return; % btnSwapLeftAndRight_Callback


% --- Executes on button press in btnDelete.
function btnDelete_Callback(hObject, eventdata, handles)
% hObject    handle to btnDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	% Get a list of those indexes that are selected, so we know which images to process.
    itemsSelected = get(handles.lstStitchedImageList, 'value');
    numberOfSelectedFiles = uint16(length(itemsSelected));
	if numberOfSelectedFiles == 0
		% Disable Delete button.
		set(handles.btnDelete, 'Enable', 'off');
		return;
	end
	% Then get list of all of the filenames in the list,
	% regardless of whether they are selected or not.
    ListOfStitchedImageNames = get(handles.lstStitchedImageList, 'string');
	promptUser = 1;
	for j=1:numberOfSelectedFiles    % Loop though all selected indexes.
        index = itemsSelected(j);    % Get the next selected index.
        % Get the filename for this selected index.
        fileName = strcat(cell2mat(ListOfStitchedImageNames(index)));
		fullFileName = fullfile(handles.stitchedImageFolder, fileName);
		if exist(fullFileName, 'file')
			% If it exists, delete it, after prompting user.
			promptString = sprintf('Delete file:\n %s?', fullFileName);
			button = questdlg(promptString,'Delete?', 'Yes', 'No', 'Cancel','Yes');
			% Comment out below because MATLAB doesn't support 4 options.
% 			button = questdlg(promptString,'Delete?', 'Yes', 'No', 'Cancel', 'All','Yes');
% 			if strcmpi(button, 'All')
% 				% User selected All option, so don't prompt anymore.
% 				promptUser = 0;
% 			end
			if strcmpi(button, 'Yes') || promptUser == 0
				% The action that the delete function takes on deleted files depends upon
				% the setting of the MATLAB recycle state. If you set the recycle state to on,
				% MATLAB moves deleted files to your recycle bin or temporary directory. With
				% the recycle state set to off (the default), deleted files
				% are permanently removed from the system.
				recycle on;
				delete(fullFileName);
			elseif strcmpi(button, 'Cancel')
				break;
			end
		end		
	end
		% Update the stitched image listbox.
	handles = LoadImageList(handles, 0, 0, 1);
    guidata(hObject, handles);
	return; % btnDelete_Callback


%=====================================================================
% --- Executes on mouse press on right image.
function Right_Image_Click_Event(src, eventdata)
	global rightImage;
	if ~isempty(rightImage)
		imtool(rightImage);
	end
	return; % Right_Image_Click_Event

%=====================================================================
% --- Executes on mouse press over axes background.
function Left_Image_Click_Event(src, eventdata)
% hObject    handle to axesLeftImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	global leftImage;
	if ~isempty(leftImage)
		imtool(leftImage);
	end
	return; % Left_Image_Click_Event
	
%=====================================================================
% --- Executes on mouse press on right image.
function Stitched_Image_Click_Event(src, eventdata)
	global stitchedImage;
	if ~isempty(stitchedImage)
		imtool(stitchedImage);
	end
	return; % Stitched_Image_Click_Event



% --- Executes on button press in chkPromptForFileName.
function chkPromptForFileName_Callback(hObject, eventdata, handles)
% hObject    handle to chkPromptForFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chkPromptForFileName
	fileNameOption = get(handles.chkPromptForFileName, 'value');
	if fileNameOption == 1
		% It calls uiputfile, which always prompts for overwriting.
		set(handles.chkPromptToOverwrite, 'value', 1);
	end
    guidata(hObject, handles);
	return; % chkPromptForFileName_Callback

% --- Executes on button press in chkPromptToOverwrite.
function chkPromptToOverwrite_Callback(hObject, eventdata, handles)
% hObject    handle to chkPromptToOverwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chkPromptToOverwrite
	overwriteOption = get(handles.chkPromptToOverwrite, 'value');
	fileNameOption = get(handles.chkPromptForFileName, 'value');
	if fileNameOption == 1 && overwriteOption == 0
		% It calls uiputfile, which always prompts for overwriting.
		% So re-check this option, don't let them uncheck it.
		set(handles.chkPromptToOverwrite, 'value', 1);
		message = sprintf('As long as the Prompt for file name option is checked\nit will always prompt before overwriting because\nthat is forced by the MATLAB uiputfile() function.\nTo uncheck this, you must uncheck the Prompt for file name option first.');
		msgboxw(message);
	end
    guidata(hObject, handles);
	return; % chkPromptToOverwrite_Callback


% --- Executes on selection change in popMode.
function popMode_Callback(hObject, eventdata, handles)
% hObject    handle to popMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popMode
	handles = UpdateAnalyzeButtonCaption(handles);
	SetSizeMismatchPopupEntries(handles);
	return;


% --- Executes during object creation, after setting all properties.
function popMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnMakeFoldersMatch.
function btnMakeFoldersMatch_Callback(hObject, eventdata, handles)
% hObject    handle to btnMakeFoldersMatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	handles.rightImageFolder = handles.leftImageFolder;
	% Load up the right image listbox only.
	handles = LoadImageList(handles, 0, 1, 0);
	guidata(hObject, handles);
	return; % from btnMakeFoldersMatch_Callback

%-------------------------------------------------------------------------------------------------------	
% Make sure none are less than 1 or greater than the canvass size.
function [topLine1 topLine2 bottomLine1 bottomLine2 leftColumn1 leftColumn2 rightColumn1 rightColumn2] = CheckBounds(expandedCanvasSize, topLine1, topLine2, bottomLine1, bottomLine2, leftColumn1, leftColumn2, rightColumn1, rightColumn2)
	if topLine1 < 1
		topLine1 = 1;
	end
	if topLine2 < 1
		topLine2 = 1;
	end
	if bottomLine1 < 1
		bottomLine1 = 1;
	end
	if bottomLine2 < 1
		bottomLine2 = 1;
	end
	if leftColumn1 < 1
		leftColumn1 = 1;
	end
	if leftColumn2 < 1
		leftColumn2 = 1;
	end
	if rightColumn1 < 1
		rightColumn1 = 1;
	end
	if rightColumn2 < 1
		rightColumn2 = 1;
	end
	
	if topLine1 > expandedCanvasSize(1)
		topLine1 = expandedCanvasSize(1);
	end
	if topLine2 > expandedCanvasSize(1)
		topLine2 = expandedCanvasSize(1);
	end
	if bottomLine1 > expandedCanvasSize(1)
		bottomLine1 = expandedCanvasSize(1);
	end
	if bottomLine2 > expandedCanvasSize(1)
		bottomLine2 = expandedCanvasSize(1);
	end
	if leftColumn1 > expandedCanvasSize(2)
		leftColumn1 = expandedCanvasSize(2);
	end
	if leftColumn2 > expandedCanvasSize(2)
		leftColumn2 = expandedCanvasSize(2);
	end
	if rightColumn1 > expandedCanvasSize(2)
		rightColumn1 = expandedCanvasSize(2);
	end
	if rightColumn2 > expandedCanvasSize(2)
		rightColumn2 = expandedCanvasSize(2);
	end
	
%-------------------------------------------------------------------------------------------------------	
% Make sure none are less than 1 because MATLAB arrays start at index 1.
function [topLine1 topLine2 leftColumn1 leftColumn2] = ClipToOne(topLine1, topLine2, leftColumn1, leftColumn2)
	if topLine1 < 1
		topLine1 = 1;
	end
	if topLine2 < 1
		topLine2 = 1;
	end
	if leftColumn1 < 1
		leftColumn1 = 1;
	end
	if leftColumn2 < 1
		leftColumn2 = 1;
	end

%=====================================================================
% Pops up a message box and waits for the user to click OK.
function msgboxw(in_strMessage)
    uiwait(msgbox(in_strMessage));
    return
