function varargout = RealspaceLattice01(varargin)
%% This is used to identify the location of Bragg peaks in the mean CBED image
%
% Colin Ophus
% National Center for Electron Microscopy, LBNL
% 2018/19/09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RealspaceLattice01_OpeningFcn, ...
                   'gui_OutputFcn',  @RealspaceLattice01_OutputFcn, ...
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

% --- Executes just before RealspaceLattice01 is made visible.
function RealspaceLattice01_OpeningFcn(hObject, eventdata, handles, varargin)

% initialize main struct "s" if not provided
if isstruct(varargin{1})
    handles.s = varargin{1};
    if nargin == 5
        I = varargin{2};
        handles.s.imageIntMean = mean(I(:));
        I = I - handles.s.imageIntMean;
        handles.s.imageIntSD = sqrt(mean(I(:).^2));
        I = I / handles.s.imageIntSD;
        handles.s.image = I;
    end
else
    handles.s = struct();
    % get image data, normalize
    I = varargin{1};
    handles.s.imageIntMean = mean(I(:));
    I = I - handles.s.imageIntMean;
    handles.s.imageIntSD = sqrt(mean(I(:).^2));
    I = I / handles.s.imageIntSD;
    handles.s.image = I;
    % inialize all needed variables
    handles.s.imageSize = size(handles.s.image);
    xr = round(handles.s.imageSize(1)*[.4 .6]) + [1 0];
    yr = round(handles.s.imageSize(2)*[.4 .6]) + [1 0];
    handles.s.p = [xr(1) yr(1);xr(2) yr(1);xr(2) yr(2);xr(1) yr(2)];
    handles.s.pos = [-1 -1 0 0 0 0];
    % Control variables
    handles.s.peakSmoothing = 1;
    handles.s.peakThresh = 0;
    handles.s.peakMinDist = 0;
    % lattice variables
    handles.s.lat = ...
        [handles.s.imageSize/2;
        [handles.s.imageSize(1)/50 0];
        [0 handles.s.imageSize(2)/50]];
    handles.s.latNumPlot = 4;
    handles.s.latSiteFrac = [1 1];
    % Strain variables
    handles.s.strainIntRange = [-1 1]; 
    handles.s.strainDispRange = [-1 1]; 
    handles.s.strainRange = [-5 5]; 
    handles.s.strainSmooth = 1;
    % Mean UC
    handles.s.UCsize = [32 32];
    handles.s.UCrep = [2 2];
    handles.s.UCsigma = 1;
    handles.s.UCmean = zeros(handles.s.UCsize);
    handles.s.UCsig = zeros(handles.s.UCsize);
    handles.s.UCcount = zeros(handles.s.UCsize);
end

if ~isfield(handles.s,'imagePlotScale')
    handles.s.imagePlotScale = [-1 3];
end

% GUI initialization
tic
handles.flagROI = 0;
handles.flagLat = 0;

% Plot x,y coordinate vectors in lower right
set(handles.figure1,'currentaxes',handles.axesXYcoord);
cla
% hold on
line([1 0 0],[0 0 1],'linewidth',2,'color','k')
line([-.2 0 -.2]+1.025,[-.2 0 .2],'linewidth',2,'color','k')
line([-.2 0 .2],[-.2 0 -.2]+1.025,'linewidth',2,'color','k')
text(1,0.35,'y','fontsize',12,'horizontalalign','center')
text(0.3,1,'x','fontsize',12)
% hold off
axis equal off
xlim([-.3 1.3])
ylim([-.3 1.3])
set(gca,'ydir','reverse')


% Image coordinate system
[handles.ya,handles.xa] = meshgrid(1:handles.s.imageSize(2),...
    1:handles.s.imageSize(1));

% Initial plotting for main axes
set(handles.figure1,'currentaxes',handles.axesMain);
imagesc(handles.s.image);
caxis(handles.s.imagePlotScale)
hold on
% Image and lattice peak positions
handles.posImage = scatter(handles.s.pos(:,2),handles.s.pos(:,1),...
    'marker','.','markeredgecolor',[1 0 0],'sizedata',80);
handles.posFit = scatter(handles.s.pos(:,6),handles.s.pos(:,5),...
    'marker','+','sizedata',80,...
    'markeredgecolor',[0 0 0],'markerfacecolor','none');
% Draw ROI selection region
handles = drawROI(handles);
% Draw lattice
handles.u = plot(handles.s.lat(1,2)...
    + [0 handles.s.lat(2,2)*handles.s.latNumPlot],...
    handles.s.lat(1,1)...
    + [0 handles.s.lat(2,1)*handles.s.latNumPlot],...
    'linewidth',2,'color','r',...
    'marker','o','markerfacecolor','w','markersize',8);
handles.v = plot(handles.s.lat(1,2)...
    + [0 handles.s.lat(3,2)*handles.s.latNumPlot],...
    handles.s.lat(1,1)...
    + [0 handles.s.lat(3,1)*handles.s.latNumPlot],...
    'linewidth',2,'color','b',...
    'marker','o','markerfacecolor','w','markersize',8);
handles.or = plot(handles.s.lat(1,2),...
    handles.s.lat(1,1),...
    'linewidth',2,'color','k',...
    'marker','o','markerfacecolor','w','markersize',10);



hold off
axis equal off
colormap(gray(256))
% handles.s.imageSize
ylim([1 handles.s.imageSize(1)])
xlim([1 handles.s.imageSize(2)])



% Set up GUI text values and drawing things
updatePeakValues(handles)
updateLatValues(handles);
updateStrainValues(handles);


% If strain maps exist, activate buttons
if isfield(handles.s,'strainEuu')
    set(handles.pushbuttonStrainExport,'enable','on');
    set(handles.pushbuttonStrainU,'enable','on');
    set(handles.pushbuttonStrainV,'enable','on');
    set(handles.pushbuttonStrainRadial,'enable','on');
    set(handles.pushbuttonStrainTheta,'enable','on');
    set(handles.pushbuttonStrainEuu,'enable','on');
    set(handles.pushbuttonStrainEvv,'enable','on');
    set(handles.pushbuttonStrainEuv,'enable','on');
    set(handles.pushbuttonQuiverPlot,'enable','on');
    set(handles.pushbuttonQuiverColour,'enable','on');
end


% intialize mean UC to black
drawMeanUC(handles)
% Update UC size values
updateUCValues(handles)

% Choose default command line output for RealspaceLattice01
handles.output = handles.s;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes RealspaceLattice01 wait for user response (see UIRESUME)
uiwait(handles.figure1);


function drawMeanUC(handles)
% mean UC plotting
set(handles.figure1,'currentaxes',handles.axesMeanUC);
imagesc(repmat(handles.s.UCmean,handles.s.UCrep))
axis equal off
colormap(gray(256))
% caxis([0 1])


% --- Outputs from this function are returned to the command line.
function varargout = RealspaceLattice01_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.s;
delete(handles.figure1);


function drawpeakStrain(handles)
% mean UC plotting
set(handles.figure1,'currentaxes',handles.axespeakStrain);
imagesc(repmat(handles.s.UCmean,handles.s.UCrep))
axis equal off
colormap(gray(256))
% caxis([0 1])


function updatePeakValues(handles)
set(handles.editPeakImageSmoothing,'string',handles.s.peakSmoothing);
set(handles.editPeakThresh,'string',handles.s.peakThresh);
set(handles.editPeakMinDist,'string',handles.s.peakMinDist);

function updateLatValues(handles)
set(handles.editLatOrX,'string',sprintf('%.02f',handles.s.lat(1,1)));
set(handles.editLatOrY,'string',sprintf('%.02f',handles.s.lat(1,2)));
set(handles.editLatUX,'string',sprintf('%.04f',handles.s.lat(2,1)));
set(handles.editLatUY,'string',sprintf('%.04f',handles.s.lat(2,2)));
set(handles.editLatVX,'string',sprintf('%.04f',handles.s.lat(3,1)));
set(handles.editLatVY,'string',sprintf('%.04f',handles.s.lat(3,2)));
% Lower
set(handles.editLatNumPlot,'string',handles.s.latNumPlot);
set(handles.editLatSiteFracU,'string',handles.s.latSiteFrac(1));
set(handles.editLatSiteFracV,'string',handles.s.latSiteFrac(2));
% Drawing
set(handles.or,'xdata',handles.s.lat(1,2));
set(handles.or,'ydata',handles.s.lat(1,1));
set(handles.u,'xdata',handles.s.lat(1,2) ...
    +[0 handles.s.lat(2,2)]*handles.s.latNumPlot);
set(handles.u,'ydata',handles.s.lat(1,1) ...
    +[0 handles.s.lat(2,1)]*handles.s.latNumPlot);
set(handles.v,'xdata',handles.s.lat(1,2) ...
    +[0 handles.s.lat(3,2)]*handles.s.latNumPlot);
set(handles.v,'ydata',handles.s.lat(1,1) ...
    +[0 handles.s.lat(3,1)]*handles.s.latNumPlot);
% strain lattice
uS = handles.s.lat(2,:);
uS = uS / norm(uS);
set(handles.textStrainLatUX,'string',sprintf('%.04f',uS(1)));
set(handles.textStrainLatUY,'string',sprintf('%.04f',uS(2)));
set(handles.textStrainLatVX,'string',sprintf('%.04f',-uS(2)));
set(handles.textStrainLatVY,'string',sprintf('%.04f',uS(1)));



function updateStrainValues(handles)
set(handles.editStrainIntMin,'string',handles.s.strainIntRange(1));
set(handles.editStrainIntMax,'string',handles.s.strainIntRange(2));
set(handles.editStrainDispMin,'string',handles.s.strainDispRange(1));
set(handles.editStrainDispMax,'string',handles.s.strainDispRange(2));
set(handles.editStrainMin,'string',handles.s.strainRange(1));
set(handles.editStrainMax,'string',handles.s.strainRange(2));
set(handles.editStrainSmooth,'string',handles.s.strainSmooth);






function handles = drawROI(handles)
if isfield(handles,'roi')
    delete(handles.roi);
end
handles.roi = patch(handles.s.p(:,2),handles.s.p(:,1),'r',...
    'linewidth',2,'edgecolor','g','facecolor','none',...
    'marker','s','markerfacecolor','w','markersize',8);


% function handles = drawLat(handles)
% if isfield(handles,'or')
%     delete(handles.or);
% end
% if isfield(handles,'u')
%     delete(handles.u);
% end
% if isfield(handles,'v')
%     delete(handles.v);
% end



% set(handles.roi,'vertices',[handles.s.p(:,2) handles.s.p(:,1)]);
% Update handles structure
% guidata(hObject, handles);




function editPeakMinDist_Callback(hObject, eventdata, handles)
% Set peak image smoothing value
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='.') == 1
    t = str2double(t);  
    handles.s.peakMinDist = t;
end
updatePeakValues(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editPeakMinDist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPeakThresh_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.peakThresh = t;
end
updatePeakValues(handles);
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function editPeakThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPeakImageSmoothing_Callback(hObject, eventdata, handles)
% Set peak image smoothing value
t = get(hObject,'String');
% if min(isstrprop(t,'digit')) == 1
if min(isstrprop(t,'digit') | t=='.') == 1
    t = min(str2double(t),16);  % Force max to 16
    handles.s.peakSmoothing = t;
end
updatePeakValues(handles);
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function editPeakImageSmoothing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPeakImageSmoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPeakMax.
function pushbuttonPeakMax_Callback(hObject, eventdata, handles)
% Use LP filtering
[qx,qy] = makeFourierCoords(handles.s.imageSize,1/norm(handles.s.lat(2,:)));
[qya,qxa] = meshgrid(qy,qx);
q2 = qxa.^2 + qya.^2;
wFilt = 1-1./sqrt(1+q2.^8/(.5^16));  % Order 8 Butterworth filter
I = handles.s.image;
I = I - mean(I(:));
I = real(ifft2(fft2(I).*wFilt));
% Scale image by SDs in local area
in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
I = I - mean(I(in));
I = I / sqrt(mean(I(in).^2));
% find peaks from image maxima
if handles.s.peakSmoothing > 0
    sm = fspecial('gaussian',2*ceil(3*handles.s.peakSmoothing),...
        handles.s.peakSmoothing);
    Ism = conv2(I,sm,'same');
else
    Ism = I;
end
p =   Ism > circshift(Ism,[-1 -1]) ...
    & Ism > circshift(Ism,[ 0 -1]) ...
    & Ism > circshift(Ism,[ 1 -1]) ...
    & Ism > circshift(Ism,[-1  0]) ...
    & Ism > circshift(Ism,[ 1  0]) ...
    & Ism > circshift(Ism,[-1  1]) ...
    & Ism > circshift(Ism,[ 0  1]) ...
    & Ism > circshift(Ism,[ 1  1]) ...
    & Ism > handles.s.peakThresh;
[xp,yp,Ip] = find(p.*Ism);
% If min distance is > 0, filter out peaks
if handles.s.peakMinDist > 0
    data = sortrows([xp yp Ip],3);
    del = false(size(data,1),1);
    r2 = handles.s.peakMinDist^2;
    for a0 = 1:(size(data,1)-1)
        d2 = (data(a0,1)-data((a0+1):end,1)).^2 ...
            + (data(a0,2)-data((a0+1):end,2)).^2;
        if min(d2) < r2
            del(a0) = true;
        end
    end
    data(del,:) = [];  
    xp = data(:,1);
    yp = data(:,2);
end
% Filter via ROI
in = inpolygon(xp,yp,handles.s.p(:,1),handles.s.p(:,2));
xp(~in) = [];
yp(~in) = [];
% Output values
handles.s.pos = [xp yp zeros(length(xp),4)];
set(handles.posImage,'xdata',handles.s.pos(:,2),...
    'ydata',handles.s.pos(:,1));
set(handles.posFit,'xdata',-1,'ydata',-1);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbuttonPeakMin.
function pushbuttonPeakMin_Callback(hObject, ~, handles)
% Use LP filtering
[qx,qy] = makeFourierCoords(handles.s.imageSize,1/norm(handles.s.lat(2,:)));
[qya,qxa] = meshgrid(qy,qx);
q2 = qxa.^2 + qya.^2;
wFilt = 1-1./sqrt(1+q2.^8/(.5^16));  % Order 8 Butterworth filter
I = handles.s.image;
I = I - mean(I(:));
I = real(ifft2(fft2(I).*wFilt));
% % Scale image by SDs in local area
in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
I = I - mean(I(in));
I = I / sqrt(mean(I(in).^2));
% find peaks from image minima
if handles.s.peakSmoothing > 0
    sm = fspecial('gaussian',2*ceil(3*handles.s.peakSmoothing),...
        handles.s.peakSmoothing);
    Ism = -conv2(I,sm,'same');
else
    Ism = -I;
end
% figure(1)
% clf
% imagesc(Ism)
% axis equal off
% colormap(violetFire)
p =   Ism > circshift(Ism,[-1 -1]) ...
    & Ism > circshift(Ism,[ 0 -1]) ...
    & Ism > circshift(Ism,[ 1 -1]) ...
    & Ism > circshift(Ism,[-1  0]) ...
    & Ism > circshift(Ism,[ 1  0]) ...
    & Ism > circshift(Ism,[-1  1]) ...
    & Ism > circshift(Ism,[ 0  1]) ...
    & Ism > circshift(Ism,[ 1  1]) ...
    & Ism > handles.s.peakThresh;
[xp,yp,Ip] = find(p.*Ism);
% If min distance is > 0, filter out peaks
if handles.s.peakMinDist > 0
    data = sortrows([xp yp Ip],3);
    del = false(size(data,1),1);
    r2 = handles.s.peakMinDist^2;
    for a0 = 1:(size(data,1)-1)
        d2 = (data(a0,1)-data((a0+1):end,1)).^2 ...
            + (data(a0,2)-data((a0+1):end,2)).^2;
        if min(d2) < r2
            del(a0) = true;
        end
    end
    data(del,:) = [];  
    xp = data(:,1);
    yp = data(:,2);
end
% Filter via ROI
in = inpolygon(xp,yp,handles.s.p(:,1),handles.s.p(:,2));
xp(~in) = [];
yp(~in) = [];
% Output values
handles.s.pos = [xp yp zeros(length(xp),4)];
set(handles.posImage,'xdata',handles.s.pos(:,2),...
    'ydata',handles.s.pos(:,1));
set(handles.posFit,'xdata',-1,'ydata',-1);
% Update handles structure
guidata(hObject, handles);





function editLatOrX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatOrX_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    xR = [max(2-[0 handles.s.lat(2,1) handles.s.lat(3,1)] ...
        *handles.s.latNumPlot) ...
        min(handles.s.imageSize(1)-1 ...
        -[0 handles.s.lat(2,1) handles.s.lat(3,1)] ...
        *handles.s.latNumPlot)];
    t = max(min(t,xR(2)),xR(1));
    handles.s.lat(1,1) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function editLatOrY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatOrY_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    yR = [max(2-[0 handles.s.lat(2,2) handles.s.lat(3,2)] ...
        *handles.s.latNumPlot) ...
        min(handles.s.imageSize(2)-1 ...
        -[0 handles.s.lat(2,2) handles.s.lat(3,2)] ...
        *handles.s.latNumPlot)];
    t = max(min(t,yR(2)),yR(1));
    handles.s.lat(1,2) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function editLatUX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatUX_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.lat(2,1) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);

function editLatUY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatUY_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.lat(2,2) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);

function editLatVX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatVX_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.lat(3,1) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);

function editLatVY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatVY_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.lat(3,2) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function editLatNumPlot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatNumPlot_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),1);
    handles.s.latNumPlot = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function editLatSiteFracU_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatSiteFracU_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),1);
    handles.s.latSiteFrac(1) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function editLatSiteFracV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editLatSiteFracV_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),1);
    handles.s.latSiteFrac(2) = t;
end
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);


function pushbuttonLatRefine_Callback(hObject, eventdata, handles)
% Need at least 3 points to fit linear 2D lattice
if size(handles.s.pos,1) >= 3
    f = handles.s.latSiteFrac;
    or = handles.s.lat(1,:);
    u = handles.s.lat(2,:);
    v = handles.s.lat(3,:);
    
    dMax = sqrt(max((handles.s.pos(:,1)-or(1)).^2 + (handles.s.pos(:,2)-or(2)).^2));
    dMin = sqrt(min((handles.s.pos(:,1)-or(1)).^2 ...
        + (handles.s.pos(:,2)-or(2)).^2)) + 5*mean([norm(u) norm(v)]);
    
    for r = dMin:20:dMax
        sub = (handles.s.pos(:,1) - or(1)).^2 ...
            + (handles.s.pos(:,2) - or(2)).^2 < r^2;
        
        % Compute a,b values
        a = round(f(1)*((handles.s.pos(sub,2)-or(2))*v(1) ...
            - (handles.s.pos(sub,1)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))/f(1);
        b = round(f(2)*((handles.s.pos(sub,2)-or(2))*u(1) ...
            - (handles.s.pos(sub,1)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))/f(2);
        % Refine lattice
        A = [ones(sum(sub),1) a b];
        xbeta = A \ handles.s.pos(sub,1);
        ybeta = A \ handles.s.pos(sub,2);
        or = [xbeta(1) ybeta(1)];
        u = [xbeta(2) ybeta(2)];
        v = [xbeta(3) ybeta(3)];
    end
    
    % Compute a,b values
    a = round(f(1)*((handles.s.pos(:,2)-or(2))*v(1) ...
        - (handles.s.pos(:,1)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))/f(1);
    b = round(f(2)*((handles.s.pos(:,2)-or(2))*u(1) ...
        - (handles.s.pos(:,1)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))/f(2);
    % Refine lattice
    A = [ones(size(handles.s.pos,1),1) a b];
    xbeta = A \ handles.s.pos(:,1);
    ybeta = A \ handles.s.pos(:,2);
    or = [xbeta(1) ybeta(1)];
    u = [xbeta(2) ybeta(2)];
    v = [xbeta(3) ybeta(3)];
    
    handles.s.lat(1,:) = or;
    handles.s.lat(2,:) = u;
    handles.s.lat(3,:) = v;
    % generated fitted positions
    xf = or(1) + a*u(1) + b*v(1);
    yf = or(2) + a*u(2) + b*v(2);
    % Export and plot fitted positions
    handles.s.pos(:,3:6) = [a b xf yf];
    set(handles.posFit,'xdata',yf,'ydata',xf);
    updateLatValues(handles);
    % Update handles structure
    guidata(hObject, handles);
end





% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% handles.output = handles.s;
% delete(hObject);





% ******** GUI CONTROL FUNCTIONS **********
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% Get cursor position
xyMain = get(handles.axesMain,'Currentpoint');
xy = [xyMain(1,2) xyMain(1,1)];
% Check to see if xy is on top of a current roi point
minDist = (.015*handles.s.imageSize(1));
[val,ind] = min((handles.s.p(:,1)-xy(1)).^2 ...
    + (handles.s.p(:,2)-xy(2)).^2);
% Determine whether click is single or double click
if toc < 1/3  
    % double click detected
    if val < minDist^2
        % point deletion if possible (more than 3 points)
        if size(handles.s.p,1) > 3
            handles.s.p(ind,:) = [];
        end
        handles = drawROI(handles); 
    else
        % point creation
        testRTuv = zeros(size(handles.s.p,1),6);
        for a0 = 1:size(handles.s.p,1)
            u = handles.s.p(mod(a0,size(handles.s.p,1))+1,:) ...
                - handles.s.p(a0,:);
            v = [u(2) u(1)] / norm(u);
            % Solve click location as superpositions of u and v
            T = ((xy(2)-handles.s.p(a0,2))*v(1) ...
                - (xy(1)-handles.s.p(a0,1))*v(2))/(v(1)*u(2)-v(2)*u(1));
            R = ((xy(2)-handles.s.p(a0,2))*u(1) ...
                - (xy(1)-handles.s.p(a0,1))*u(2))/(v(2)*u(1)-v(1)*u(2));
            testRTuv(a0,:) = [T R u v];
        end
        [minR,ind] = min(abs(testRTuv(:,2)));
        if minR < minDist && testRTuv(ind,1) > 0 && testRTuv(ind,1) < 1
            % Generate new point
            pNew = handles.s.p(ind,:) ...
                + testRTuv(ind,1)*testRTuv(ind,3:4);
            if ind == size(handles.s.p,1)
                pp = [handles.s.p; pNew];
            else
                pp = [handles.s.p(1:ind,:);
                    pNew; handles.s.p((ind+1):end,:)];
            end
            handles.s.p = pp;
            handles = drawROI(handles); 
        end
    end
else
    % single click case --> move points or lattice vectors
    % Check distance to  ROI points
    if val < minDist^2
        handles.flagROI = ind;
    else    
        % else check for distance to lattice points
        lat = handles.s.lat;
        lat(2,:) = lat(1,:) + lat(2,:)*handles.s.latNumPlot;
        lat(3,:) = lat(1,:) + lat(3,:)*handles.s.latNumPlot;
        [val,ind] = min((lat(:,1)-xy(1)).^2 ...
            + (lat(:,2)-xy(2)).^2);
        if val < minDist^2
            handles.flagLat = ind;
        end
    end
    
end
tic
% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
if handles.flagROI > 0
    xyMain = get(handles.axesMain,'Currentpoint');
    xy = [xyMain(1,2) xyMain(1,1)];
    if xy(1) < 1
        xy(1) = 1;
    end
    if xy(2) < 1
        xy(2) = 1;
    end
    if xy(1) > handles.s.imageSize(1)
        xy(1) = handles.s.imageSize(1);
    end
    if xy(2) > handles.s.imageSize(2)
        xy(2) = handles.s.imageSize(2);
    end
    handles.s.p(handles.flagROI,:) = xy;
    handles = drawROI(handles);
end
if handles.flagLat > 0
    xyMain = get(handles.axesMain,'Currentpoint');
    xy = [xyMain(1,2) xyMain(1,1)];
    if handles.flagLat == 1
        xR = [max(2-[0 handles.s.lat(2,1) handles.s.lat(3,1)] ...
            *handles.s.latNumPlot) ...
            min(handles.s.imageSize(1)-1 ...
            -[0 handles.s.lat(2,1) handles.s.lat(3,1)] ...
            *handles.s.latNumPlot)];
        yR = [max(2-[0 handles.s.lat(2,2) handles.s.lat(3,2)] ...
            *handles.s.latNumPlot) ...
            min(handles.s.imageSize(2)-1 ...
            -[0 handles.s.lat(2,2) handles.s.lat(3,2)] ...
            *handles.s.latNumPlot)];
    else
        xR = [1 handles.s.imageSize(1)];
        yR = [1 handles.s.imageSize(2)];
    end
    if xy(1) < xR(1)
        xy(1) = xR(1);
    end
    if xy(2) < yR(1)
        xy(2) = yR(1);
    end
    if xy(1) > xR(2)
        xy(1) = xR(2);
    end
    if xy(2) > yR(2)
        xy(2) = yR(2);
    end
    if handles.flagLat == 1
        handles.s.lat(1,:) = xy;
        
    elseif handles.flagLat == 2
        handles.s.lat(2,:) = (xy - handles.s.lat(1,:)) ...
            / handles.s.latNumPlot;
    else
        handles.s.lat(3,:) = (xy - handles.s.lat(1,:)) ...
            / handles.s.latNumPlot;
    end
    updateLatValues(handles);
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
if handles.flagROI > 0
    handles.flagROI = 0;
end
if handles.flagLat > 0
    handles.flagLat = 0;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbuttonImageZoom.
function pushbuttonImageZoom_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonImageZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbuttonImageZoom
if get(hObject,'Value') == true
    set(handles.figure1,'currentaxes',handles.axesMain);
    zoom on
else
    set(handles.figure1,'currentaxes',handles.axesMain);
    zoom off
end

function pushbuttonImageReset_Callback(hObject, eventdata, handles)
set(handles.figure1,'currentaxes',handles.axesMain);
ylim([1 handles.s.imageSize(1)])
xlim([1 handles.s.imageSize(2)])


function pushbuttonLatDynamic_Callback(hObject, eventdata, handles)


if size(handles.s.pos,1) >= 3
    f = handles.s.latSiteFrac;
    or = handles.s.lat(1,:);
    u = handles.s.lat(2,:);
    v = handles.s.lat(3,:);
    pos = handles.s.pos;
    
    radiusThresh2 = (min([norm(u) norm(v)]./f)*.25) ^ 2;
    abVals = [ones(size(pos,1),2)*-1e8 zeros(size(pos,1),1)];
    % Get origin
    [~,ind] = min((pos(:,1)-or(1)).^2 + (pos(:,2)-or(2)).^2);
    abVals(ind,:) = [0 0 1];
    % List of points not yet assigned - [ind,x,y]
    xyList = [(1:size(pos,1))' pos];  
    xyList(ind,:) = [];
    % List of trial points - [a,b,x,y]
    if min(handles.s.latSiteFrac) == 1
        dab = [-1 0;0 -1;1 0;0 1];
    else
        dab = [-1 0;0 -1;1 0;0 1;
            -1 -1;1 -1;-1 1;1 1];
        dab(:,1) = dab(:,1) / handles.s.latSiteFrac(1);
        dab(:,2) = dab(:,2) / handles.s.latSiteFrac(2);
    end
    dxy = dab(:,1)*u + dab(:,2)*v;
    xyTrial = [dab repmat(pos(ind,1:2),[size(dab,1) 1])+dxy];
    
    % While there are still trial points, keep checking!
    while size(xyTrial,1) > 0
        
        xyTrialNew = [];
        
        for a0 = 1:size(xyTrial,1)
            [val,ind] = min((xyList(:,2)-xyTrial(a0,3)).^2 ...
                + (xyList(:,3)-xyTrial(a0,4)).^2);
            
            if val < radiusThresh2
                % Write out (a,b) value
                indSelect = xyList(ind,1);
                abVals(indSelect,:) = [xyTrial(a0,1:2) 1];
                
                % Add trial points
                if isempty(xyTrialNew)
                    xyTrialNew = [repmat(xyTrial(a0,1:2),[size(dab,1) 1]) + dab ...
                        repmat(xyList(ind,2:3),[size(dab,1) 1]) + dxy];
                else
                    xyTrialNew = [xyTrialNew;
                        repmat(xyTrial(a0,1:2),[size(dab,1) 1]) + dab ...
                        repmat(xyList(ind,2:3),[size(dab,1) 1]) + dxy];
                end
                
                % Delete point from list
                xyList(ind,:) = [];
            end
        end
        
        % Remove points already found
        if ~isempty(xyTrialNew)
            [~,~,inds] = intersect(abVals(:,1:2),xyTrialNew(:,1:2),'rows');
            xyTrialNew(inds,:) = [];
        end
        
        % Merge trial points
        if ~isempty(xyTrialNew)
            [abUnique,~,inds] = unique(xyTrialNew(:,1:2),'rows');
            xyTrialNew2 = zeros(size(abUnique,1),4);
            for a0 = 1:size(abUnique,1)
                xyTrialNew2(a0,:) = [abUnique(a0,:) ...
                    mean(xyTrialNew(inds==a0,3:4),1)];
            end
        else
            xyTrialNew2 = [];
        end
        
        % Remove any a,b values that have already been detected from
        % xyTrial - protection against double-identities!
        
        
        xyTrial = xyTrialNew2;
    end

    % Delete positions not found
    keep = abVals(:,3) == 1;
    handles.s.pos = handles.s.pos(keep,:);
    a = abVals(keep,1);
    b = abVals(keep,2);
    
    % Compute a,b values
    %     a = round(f(1)*((handles.s.pos(:,2)-or(2))*v(1) ...
    %         - (handles.s.pos(:,1)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))/f(1);
    %     b = round(f(2)*((handles.s.pos(:,2)-or(2))*u(1) ...
    %         - (handles.s.pos(:,1)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))/f(2);
    % Refine lattice
    A = [ones(size(handles.s.pos,1),1) a b];
    xbeta = A \ handles.s.pos(:,1);
    ybeta = A \ handles.s.pos(:,2);
    or = [xbeta(1) ybeta(1)];
    u = [xbeta(2) ybeta(2)];
    v = [xbeta(3) ybeta(3)];
    
    % generated fitted positions
    xf = or(1) + a*u(1) + b*v(1);
    yf = or(2) + a*u(2) + b*v(2);
    % Export and plot fitted positions
    handles.s.pos(:,3:6) = [a b xf yf];
    set(handles.posFit,'xdata',yf,'ydata',xf);
    updateLatValues(handles);
    
    % Update lattice
    handles.s.lat(1,:) = or;
    handles.s.lat(2,:) = u;
    handles.s.lat(3,:) = v;
    
    % Update handles structure
    guidata(hObject, handles);
end



function editStrainIntMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainIntMin_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainIntRange(1) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainIntMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainIntMax_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainIntRange(2) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainDispMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainDispMin_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainDispRange(1) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainDispMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainDispMax_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainDispRange(2) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainMin_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainRange(1) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainMax_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='-' | t=='.') == 1
    t = str2double(t);
    handles.s.strainRange(2) = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);


function editStrainSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editStrainSmooth_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='.') == 1
    t = str2double(t);
    handles.s.strainSmooth = t;
end
updateStrainValues(handles);
% Update handles structure
guidata(hObject, handles);



function pushbuttonStrainCompute_Callback(hObject, eventdata, handles)
% Compute displacement and then strain maps
% Strain coordinate system vectors
uS = handles.s.lat(2,:);
uS = uS / norm(uS);
vS = [-uS(2) uS(1)];
%  Diplacement values
dx = handles.s.pos(:,1) - handles.s.pos(:,5);
dy = handles.s.pos(:,2) - handles.s.pos(:,6);
du = dx*uS(1) + dy*uS(2);
dv = dx*vS(1) + dy*vS(2);
%  Diplacement maps
inds = sub2ind(handles.s.imageSize,...
    round(handles.s.pos(:,1)),round(handles.s.pos(:,2)));
handles.s.strainDispU = zeros(handles.s.imageSize);
handles.s.strainDispV = zeros(handles.s.imageSize);
count = zeros(handles.s.imageSize);
handles.s.strainDispU(inds) = du;
handles.s.strainDispV(inds) = dv;
count(inds) = 1;
% apply smoothing
sigma = handles.s.strainSmooth * norm(handles.s.lat(2,:));
sm = fspecial('gaussian',2*ceil(3*sigma)+1,sigma);
handles.s.strainDispU = convolve2(handles.s.strainDispU,sm,'wrap');
handles.s.strainDispV = convolve2(handles.s.strainDispV,sm,'wrap');
count = convolve2(count,sm,'wrap');
sub = count > 0;
handles.s.strainDispU(sub) = handles.s.strainDispU(sub)./count(sub);
handles.s.strainDispV(sub) = handles.s.strainDispV(sub)./count(sub);
% Compute strain field
dUdx = (circshift(handles.s.strainDispU,[1 0]) ...
    - circshift(handles.s.strainDispU,[-1 0])) / 2;
dUdy = (circshift(handles.s.strainDispU,[0 1]) ...
    - circshift(handles.s.strainDispU,[0 -1])) / 2;
dVdx = (circshift(handles.s.strainDispV,[1 0]) ...
    - circshift(handles.s.strainDispV,[-1 0])) / 2;
dVdy = (circshift(handles.s.strainDispV,[0 1]) ...
    - circshift(handles.s.strainDispV,[0 -1])) / 2;
handles.s.strainEuu = dUdx*uS(1) + dUdy*uS(2);
handles.s.strainEvv = dVdx*vS(1) + dVdy*vS(2);
handles.s.strainEuv = (dUdx*vS(1) + dUdy*vS(2) ...
    + dVdx*uS(1) + dVdy*uS(2))/2;
% Scale strain displacement maps to be in unit cell units
handles.s.strainDispU = handles.s.strainDispU / norm(handles.s.lat(2,:));
handles.s.strainDispV = handles.s.strainDispV / norm(handles.s.lat(2,:));
% Make intensity map for scaling
handles.s.strainIval = (handles.s.image - handles.s.strainIntRange(1)) ...
    / (handles.s.strainIntRange(2) - handles.s.strainIntRange(1));
handles.s.strainIval(handles.s.strainIval<0) = 0;
handles.s.strainIval(handles.s.strainIval>1) = 1;
% Make mask for colouring
m = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
handles.s.mask = bwdist(~m) / sigma;
handles.s.mask(handles.s.mask<0) = 0;
handles.s.mask(handles.s.mask>1) = 1;
handles.s.mask = 1 - handles.s.mask;
% Enable plotting buttons
set(handles.pushbuttonStrainExport,'enable','on');
set(handles.pushbuttonStrainU,'enable','on');
set(handles.pushbuttonStrainV,'enable','on');
set(handles.pushbuttonStrainRadial,'enable','on');
set(handles.pushbuttonStrainTheta,'enable','on');
set(handles.pushbuttonStrainEuu,'enable','on');
set(handles.pushbuttonStrainEvv,'enable','on');
set(handles.pushbuttonStrainEuv,'enable','on');
set(handles.pushbuttonQuiverPlot,'enable','on');
set(handles.pushbuttonQuiverColour,'enable','on');
% Update handles structure
guidata(hObject, handles);



function pushbuttonStrainExport_Callback(hObject, eventdata, handles)
% Create and save images
fbase = inputdlg('File name stem for strain images?',...
    'Export strain map figures');
% Raw image
fname = [fbase{1} '_intensity_' ...
    num2str(handles.s.strainIntRange(1)) '_to_' ...
    num2str(handles.s.strainIntRange(2)) '_SD.png'];
imwrite(uint8(255*handles.s.strainIval),fname,'png');
% U disp
Imap = (handles.s.strainDispU - handles.s.strainDispRange(1)) ...
    / (handles.s.strainDispRange(2) - handles.s.strainDispRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_u_disp_' ...
    num2str(handles.s.strainDispRange(1)) '_to_' ...
    num2str(handles.s.strainDispRange(2)) '_UCs.png'];
imwrite(uint8(255*Istrain),fname,'png');
% V disp
Imap = (handles.s.strainDispV - handles.s.strainDispRange(1)) ...
    / (handles.s.strainDispRange(2) - handles.s.strainDispRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_v_disp_' ...
    num2str(handles.s.strainDispRange(1)) '_to_' ...
    num2str(handles.s.strainDispRange(2)) '_UCs.png'];
imwrite(uint8(255*Istrain),fname,'png');
% Radial disp
strainDispRadial = sqrt(handles.s.strainDispU.^2 ...
    +handles.s.strainDispV.^2);
Imap = (strainDispRadial) ...
    / (handles.s.strainDispRange(2));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_radial_disp_' ...
    '0_to_' num2str(handles.s.strainDispRange(2)) '_UCs.png'];
imwrite(uint8(255*Istrain),fname,'png');
% theta disp
strainTheta = atan2(handles.s.strainDispV,handles.s.strainDispU);
Imap = mod((strainTheta + 0)/(2*pi),1);
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
Istrain = zeros(size(Imap,1),size(Imap,2),3);
Istrain(:,:,1) = (1-Imap)*(2/3);
Istrain(:,:,2) = 1-handles.s.mask;
Istrain(:,:,3) = handles.s.strainIval;
Istrain = hsv2rgb(Istrain);
Istrain(Istrain<0) = 0;
Istrain(Istrain>1) = 1;
fname = [fbase{1} '_theta_disp_2_Pi_rad.png'];
imwrite(uint8(255*Istrain),fname,'png');
% Color radial + angle map
strainTheta = atan2(handles.s.strainDispV,handles.s.strainDispU);
Imap = mod((strainTheta + 0)/(2*pi),1);
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
strainDispRadial = sqrt(handles.s.strainDispU.^2 ...
    +handles.s.strainDispV.^2);
Imap2 = (strainDispRadial) ...
    / (handles.s.strainDispRange(2));
Imap2(Imap<0) = 0;
Imap2(Imap2>1) = 1;
Istrain = zeros(size(Imap,1),size(Imap,2),3);
Istrain(:,:,1) = Imap;
Istrain(:,:,2) = (1-handles.s.mask);
Istrain(:,:,3) = Imap2.*(1-handles.s.mask);
Istrain = hsv2rgb(Istrain);
Istrain(Istrain<0) = 0;
Istrain(Istrain>1) = 1;
fname = [fbase{1} '_polar_disp_' ...
    '0_to_' num2str(handles.s.strainDispRange(2)) ...
    '_UCs.png'];imwrite(uint8(255*Istrain),fname,'png');
imwrite(uint8(255*Istrain),fname,'png');


% Euu
Imap = (handles.s.strainEuu*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_e_uu_' ...
    num2str(handles.s.strainRange(1)) '_to_' ...
    num2str(handles.s.strainRange(2)) '_percent.png'];
imwrite(uint8(255*Istrain),fname,'png');
% Evv
Imap = (handles.s.strainEvv*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_e_vv_' ...
    num2str(handles.s.strainRange(1)) '_to_' ...
    num2str(handles.s.strainRange(2)) '_percent.png'];
imwrite(uint8(255*Istrain),fname,'png');
% Euv
Imap = (handles.s.strainEuv*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_e_uv_' ...
    num2str(handles.s.strainRange(1)) '_to_' ...
    num2str(handles.s.strainRange(2)) '_percent.png'];
imwrite(uint8(255*Istrain),fname,'png');
% Dilation
Idilate = ((1+handles.s.strainEuu).*(1+handles.s.strainEvv))-1;
Imap = (Idilate*100 - handles.s.strainRange(1)/2) ...
    / (handles.s.strainRange(2)/2 - handles.s.strainRange(1)/2);
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
fname = [fbase{1} '_scaling_' ...
    num2str(handles.s.strainRange(1)/2) '_to_' ...
    num2str(handles.s.strainRange(2)/2) '_percent.png'];
imwrite(uint8(255*Istrain),fname,'png');



function pushbuttonStrainU_Callback(hObject, eventdata, handles)
Imap = (handles.s.strainDispU - handles.s.strainDispRange(1)) ...
    / (handles.s.strainDispRange(2) - handles.s.strainDispRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)


function pushbuttonStrainV_Callback(hObject, eventdata, handles)
Imap = (handles.s.strainDispV - handles.s.strainDispRange(1)) ...
    / (handles.s.strainDispRange(2) - handles.s.strainDispRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)


function pushbuttonStrainRadial_Callback(hObject, eventdata, handles)
strainRadial = sqrt(handles.s.strainDispU.^2 ...
    + handles.s.strainDispV.^2);
Imap = (strainRadial) ...
    / (handles.s.strainDispRange(2));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)
function pushbuttonStrainTheta_Callback(hObject, eventdata, handles)
strainTheta = atan2(handles.s.strainDispV,handles.s.strainDispU);
Imap = mod((strainTheta + 0)/(2*pi),1);
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
% Make hue map for theta
Istrain = zeros(size(Imap,1),size(Imap,2),3);
Istrain(:,:,1) = (1-Imap)*(2/3);
Istrain(:,:,2) = 1-handles.s.mask;
Istrain(:,:,3) = handles.s.strainIval;
Istrain = hsv2rgb(Istrain);
Istrain(Istrain<0) = 0;
Istrain(Istrain>1) = 1;
plotImage(Istrain)


function pushbuttonStrainEuu_Callback(hObject, eventdata, handles)
Imap = (handles.s.strainEuu*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)

function pushbuttonStrainEvv_Callback(hObject, eventdata, handles)
Imap = (handles.s.strainEvv*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)

function pushbuttonStrainEuv_Callback(hObject, eventdata, handles)
Imap = (handles.s.strainEuv*100 - handles.s.strainRange(1)) ...
    / (handles.s.strainRange(2) - handles.s.strainRange(1));
Imap(Imap<0) = 0;
Imap(Imap>1) = 1;
% Imap = round(255*Imap) + 1;
Istrain = makeStrainImage(handles,Imap);
plotImage(Istrain)


function Istrain = makeStrainImage(handles,Imap)
Istrain = zeros(handles.s.imageSize(1),handles.s.imageSize(2),3);
% Istrain = ind2rgb(Imap,flipud(jet(256)));
Istrain(:,:,1) = (1-Imap)*(2/3);
Istrain(:,:,2) = 1;
Istrain(:,:,3) = handles.s.strainIval;
% Istrain(Istrain<0) = 0;
% Istrain(Istrain>1) = 1;
Istrain = hsv2rgb(Istrain);
% Istrain = rgb2hsv(Istrain);
Istrain = repmat(handles.s.mask.*handles.s.strainIval,[1 1 3]) ...
    + repmat(1-handles.s.mask,[1 1 3]) .* Istrain;

% Istrain = Imap;



function plotImage(Istrain)
figure(11)
clf
imagesc(Istrain)
axis equal off
% colormap(jet(256))
set(gca,'position',[0 0 1 1])


function pushbuttonImagePan_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == true
    set(handles.figure1,'currentaxes',handles.axesMain);
    pan on
else
    set(handles.figure1,'currentaxes',handles.axesMain);
    pan off
end


% --- Executes on button press in pushbuttonQuiverPlot.
function pushbuttonQuiverPlot_Callback(hObject, eventdata, handles)
% Make a quicver plot
scale = str2double(inputdlg('Scaling of vectors',...
    '[Set to "1" for unscaled]'));
figure(11)
clf
imagesc(handles.s.strainIval)
hold on

quiver(handles.s.pos(:,6),handles.s.pos(:,5),...
    (handles.s.pos(:,2)-handles.s.pos(:,6))*scale,...
    (handles.s.pos(:,1)-handles.s.pos(:,5))*scale,...
    'autoscale','off','color','r','linewidth',1)
hold off
axis equal off
colormap(gray(256))
set(gca,'position',[0 0 1 1])


function pushbuttonQuiverColour_Callback(hObject, eventdata, handles)
% Make a color quicver plot
scale = str2double(inputdlg('Scaling of vectors',...
    '[Set to "1" for unscaled]'));
figure(11)
clf
imagesc(handles.s.strainIval)
hold on
for a0 = 1:size(handles.s.pos,1)
    xy = handles.s.pos(a0,5:6);
    dxy = handles.s.pos(a0,1:2)-handles.s.pos(a0,5:6);
    c = hsv2rgb(reshape(...
        [mod(atan2(-dxy(2),-dxy(1))/(2*pi),1) 1 1],[1 1 3]));
    line([0 dxy(2)]*scale+xy(2),[0 dxy(1)]*scale+xy(1),...
        'linewidth',1,'color',c)
end
% quiver(handles.s.pos(:,6),handles.s.pos(:,5),...
%     (handles.s.pos(:,2)-handles.s.pos(:,6))*scale,...
%     (handles.s.pos(:,1)-handles.s.pos(:,5))*scale,...
%     'autoscale','off','color','r','linewidth',1)
hold off
axis equal off
colormap(gray(256))
set(gcf,'renderer','zbuffer')
set(gca,'position',[0 0 1 1])


function [qx,qy] = makeFourierCoords(N,pSize)
% This function generates image Fourier coordinates (2D)
if length(pSize) == 1
    L = N(1:2)*pSize;
elseif length(pSize) == 2
    L = N(1:2).*pSize;
end
if mod(N(1),2) == 0
    qx = circshift(((-N(1)/2):(N(1)/2-1))/L(1),[1 -N(1)/2]);
else
    qx = circshift(((-N(1)/2+.5):(N(1)/2-.5))/L(1),[1 -N(1)/2+.5]);
end
if mod(N(2),2) == 0
    qy = circshift(((-N(2)/2):(N(2)/2-1))/L(2),[1 -N(2)/2]);
else
    qy = circshift(((-N(2)/2+.5):(N(2)/2-.5))/L(2),[1 -N(2)/2+.5]);
end


function editUCu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editUCu_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),2);
    handles.s.UCsize(1) = t;
end
updateUCValues(handles);
% Update handles structure
guidata(hObject, handles);


function editUCv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editUCv_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),2);
    handles.s.UCsize(2) = t;
end
updateUCValues(handles);
% Update handles structure
guidata(hObject, handles);


function editUCrepu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editUCrepu_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),1);
    handles.s.UCrep(1) = t;
end
updateUCValues(handles);
drawMeanUC(handles);
% Update handles structure
guidata(hObject, handles);


function editUCrepv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editUCrepv_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit')) == 1
    t = max(round(str2double(t)),1);
    handles.s.UCrep(2) = t;
end
updateUCValues(handles);
drawMeanUC(handles);
% Update handles structure
guidata(hObject, handles);




function editUCsigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editUCsigma_Callback(hObject, eventdata, handles)
t = get(hObject,'String');
if min(isstrprop(t,'digit') | t=='.') == 1
    t = max(str2double(t),0);
    handles.s.UCsigma = t;
end
updateUCValues(handles);
[handles] = computeMeanUC(handles);
drawMeanUC(handles);
% Update handles structure
guidata(hObject, handles);



function pushbuttonUCsize_Callback(hObject, eventdata, handles)
% Query user about length of longer direction
% or = handles.s.lat(1,:);
u = handles.s.lat(2,:);
v = handles.s.lat(3,:);
L = inputdlg(' Length of longest UC vector in pixels? ');
if min(isstrprop(L{1},'digit')) == 1
    L = str2double(L{1});
    if norm(u) > norm(v)
        leng = round([L L*norm(v)/norm(u)]);
        leng = max(leng,[1 1]);
    else
        leng = round([L*norm(u)/norm(v) L]);
        leng = max(leng,[1 1]);
    end
    handles.s.UCsize = leng;
    %     handles = computeMeanUCsig(handles);
    handles = computeMeanUC(handles);
    drawMeanUC(handles);
    updateUCValues(handles);
    % Update handles structure
    guidata(hObject, handles);
end


function pushbuttonUCcompute_Callback(hObject, eventdata, handles)
% compute the mean UC!
% handles = computeMeanUCsig(handles);
handles = computeMeanUC(handles);
drawMeanUC(handles);
% Update handles structure
guidata(hObject, handles);



% This piece may be reused
% function [handles] = computeMeanUCsig(handles,I)
% % ROI subset
% or = handles.s.lat(1,:);
% u = handles.s.lat(2,:);
% v = handles.s.lat(3,:);
% in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
% aInd = mod(round(handles.s.UCsize(1)*((handles.ya(in)-or(2))*v(1) ...
%     - (handles.xa(in)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))-1,...
%     handles.s.UCsize(1))+1;
% bInd = mod(round(handles.s.UCsize(2)*((handles.ya(in)-or(2))*u(1) ...
%     - (handles.xa(in)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))-1,...
%     handles.s.UCsize(2))+1;
% if nargin == 1
%     handles.s.UCsig = accumarray([aInd bInd],...
%         handles.s.image(in),handles.s.UCsize);
% else
%     handles.s.UCsig = accumarray([aInd bInd],...
%         I(in),handles.s.UCsize);
% end
% handles.s.UCcount = accumarray([aInd bInd],...
%     sum(in(:)),handles.s.UCsize);
% sub = handles.s.UCcount > 0;
% handles.s.UCmean = handles.s.UCsig ./ handles.s.UCcount



% This piece may be reused
function [handles] = computeMeanUC(handles,I)
% compute signal of unit cell
% ROI subset
or = handles.s.lat(1,:);
u = handles.s.lat(2,:);
v = handles.s.lat(3,:);
in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
aInd = mod(round(handles.s.UCsize(1)*((handles.ya(in)-or(2))*v(1) ...
    - (handles.xa(in)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))-1,...
    handles.s.UCsize(1))+1;
bInd = mod(round(handles.s.UCsize(2)*((handles.ya(in)-or(2))*u(1) ...
    - (handles.xa(in)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))-1,...
    handles.s.UCsize(2))+1;
if nargin == 1
    handles.s.UCsig = accumarray([aInd bInd],...
        handles.s.image(in),handles.s.UCsize);
else
    handles.s.UCsig = accumarray([aInd bInd],...
        I(in),handles.s.UCsize);
end
handles.s.UCcount = accumarray([aInd bInd],...
    ones(length(aInd),1),handles.s.UCsize);
% Compute mean UC as UCsignal / UCcount
sm = fspecial('gaussian',2*ceil(3*handles.s.UCsigma)+1,handles.s.UCsigma);
sig = convolve2(handles.s.UCsig,sm,'wrap');
count = convolve2(handles.s.UCcount,sm,'wrap');
sub = count > 0;
handles.s.UCmean = zeros(handles.s.UCsize);
handles.s.UCmean(sub) = sig(sub) ./ count(sub);
% Compute RMS deviation of UC
imageROImeanUC = handles.s.UCsig(sub2ind(handles.s.UCsize,aInd,bInd)) ...
    ./ handles.s.UCcount(sub2ind(handles.s.UCsize,aInd,bInd));
% UCmeanSDsig = accumarray([aInd bInd],...
%     imageROImeanUC,handles.s.UCsize);
if nargin == 1
    UCmeanSDsig = accumarray([aInd bInd],...
        (handles.s.image(in)-imageROImeanUC).^2,handles.s.UCsize);
else
    UCmeanSDsig = accumarray([aInd bInd],...
        (I(in)-imageROImeanUC).^2,handles.s.UCsize);
end
UCmeanSDsig = convolve2(UCmeanSDsig,sm,'wrap');
UCmeanSDsig(sub) = UCmeanSDsig(sub) ./ count(sub);
handles.s.UCmeanSD = sqrt(UCmeanSDsig);
% handles.s.temp = UCmeanSDsig;
% temp = zeros(handles.s.imageSize);
% temp(sub2ind(handles.s.imageSize,handles.xa(in),handles.ya(in))) = imageROImeanUC;
% handles.s.temp = temp;
% SD = accumarray([aInd bInd],...
%     sum(in(:)),handles.s.UCsize);



function updateUCValues(handles)
set(handles.editUCu,'string',handles.s.UCsize(1));
set(handles.editUCv,'string',handles.s.UCsize(2));
set(handles.editUCrepu,'string',handles.s.UCrep(1));
set(handles.editUCrepv,'string',handles.s.UCrep(2));
set(handles.editUCsigma,'string',handles.s.UCsigma);


function pushbuttonUCcomputeLP_Callback(hObject, eventdata, handles)
% Use LP filtering
[qx,qy] = makeFourierCoords(handles.s.imageSize,1/norm(handles.s.lat(2,:)));
[qya,qxa] = meshgrid(qy,qx);
q2 = qxa.^2 + qya.^2;
wFilt = 1-1./sqrt(1+q2.^8/(.5^16));  % Order 8 Butterworth filter
I = handles.s.image;
I = I - mean(I(:));
I = real(ifft2(fft2(I).*wFilt));
% Scale image by SDs in local area
% in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
% I = I - mean(I(in));
% I = I / sqrt(mean(I(in).^2));
% compute the mean UC!
% handles = computeMeanUCsig(handles,I);
handles = computeMeanUC(handles);
drawMeanUC(handles);
% Update handles structure
guidata(hObject, handles);


function pushbuttonBraggLattice_Callback(hObject, eventdata, handles)
% Search for lattice in FFT space
I = handles.s.image;
% Generate filtering mask with bwdist
in = inpolygon(handles.xa,handles.ya,handles.s.p(:,1),handles.s.p(:,2));
in([1 end],:) = false;
in(:,[1 end]) = false;
mask = bwdist(~in);
mask = sin(mask*(pi/2/max(mask(:)))).^2;
% Create FFT image
Ifft = fftshift(abs(fft2(I.*mask)));
% Ifft = Ifft2.^2;
[qx,qy] = makeFourierCoords(size(Ifft),1);
qx = fftshift(qx);
qy = fftshift(qy);
% Determine plot colour range
sm = fspecial('gaussian',21,4);
Ifft = conv2(Ifft,sm,'same');
% Ifft = abs(Ifft).^0.5;
Ifft = Ifft / max(max(Ifft(2:end,2:end)));
% Ifft = Ifft / max(Ifft(:));
Ir = [0 .25];
% % Ivals = sort(Ifft2(:));
% % inds = round([.93-0.0 .99999]*length(Ivals))+[1 0];
% % Ir = Ivals(inds);
% Get input from user 
figure(11)
clf
imagesc(Ifft)
axis equal off
colormap(hot(256))
set(gca,'position',[0 0 1 1])
% Ir
caxis(Ir)
[ym,xm] = ginput(2);
% Find largest peak in Ifft close by
r2 = 8^2;
p =   Ifft > circshift(Ifft,[-1 -1]) ...
    & Ifft > circshift(Ifft,[ 0 -1]) ...
    & Ifft > circshift(Ifft,[ 1 -1]) ...
    & Ifft > circshift(Ifft,[-1  0]) ...
    & Ifft > circshift(Ifft,[ 1  0]) ...
    & Ifft > circshift(Ifft,[-1  1]) ...
    & Ifft > circshift(Ifft,[ 0  1]) ...
    & Ifft > circshift(Ifft,[ 1  1]);
[xp,yp,Ip] = find(p.*Ifft);
sub = (xp-xm(1)).^2 + (yp-ym(1)).^2 < r2;
data = sortrows([xp(sub) yp(sub) Ip(sub)],-3);
xy = data(1,1:2);
Icut = Ifft(xy(1)+(-1:1),xy(2)+(-1:1));
dx = (Icut(3,2)-Icut(1,2))/(2*Icut(2,2)-Icut(3,2)-Icut(1,2))/2;
dy = (Icut(2,3)-Icut(2,1))/(2*Icut(2,2)-Icut(2,3)-Icut(2,1))/2;
if dx > 0
    pqx = qx(xy(1))*(1-dx) + qx(xy(1)+1)*dx;
else
    pqx = qx(xy(1))*(1+dx) + qx(xy(1)-1)*-dx;
end
if dy > 0
    pqy = qy(xy(2))*(1-dy) + qy(xy(2)+1)*dy;
else
    pqy = qy(xy(2))*(1+dy) + qy(xy(2)-1)*-dy;
end
qvec = [pqx pqy];
handles.s.lat(2,:) = qvec./norm(qvec)^2;
% uNew = qvec./norm(qvec)^2;
sub = (xp-xm(2)).^2 + (yp-ym(2)).^2 < r2;
data = sortrows([xp(sub) yp(sub) Ip(sub)],-3);
xy = data(1,1:2);
Icut = Ifft(xy(1)+(-1:1),xy(2)+(-1:1));
dx = (Icut(3,2)-Icut(1,2))/(2*Icut(2,2)-Icut(3,2)-Icut(1,2))/2;
dy = (Icut(2,3)-Icut(2,1))/(2*Icut(2,2)-Icut(2,3)-Icut(2,1))/2;
if dx > 0
    pqx = qx(xy(1))*(1-dx) + qx(xy(1)+1)*dx;
else
    pqx = qx(xy(1))*(1+dx) + qx(xy(1)-1)*-dx;
end
if dy > 0
    pqy = qy(xy(2))*(1-dy) + qy(xy(2)+1)*dy;
else
    pqy = qy(xy(2))*(1+dy) + qy(xy(2)-1)*-dy;
end
qvec = [pqx pqy];
handles.s.lat(3,:) = qvec./norm(qvec)^2;
% generated fitted positions
or = handles.s.lat(1,:);
u = handles.s.lat(2,:);
v = handles.s.lat(3,:);
% Compute a,b values
f = handles.s.latSiteFrac;
a = round(f(1)*((handles.s.pos(:,2)-or(2))*v(1) ...
    - (handles.s.pos(:,1)-or(1))*v(2))/(v(1)*u(2)-v(2)*u(1)))/f(1);
b = round(f(2)*((handles.s.pos(:,2)-or(2))*u(1) ...
    - (handles.s.pos(:,1)-or(1))*u(2))/(v(2)*u(1)-v(1)*u(2)))/f(2);
xf = or(1) + a*u(1) + b*v(1);
yf = or(2) + a*u(2) + b*v(2);
% Export and plot fitted positions
handles.s.pos(:,3:6) = [a b xf yf];
set(handles.posFit,'xdata',yf,'ydata',xf);
updateLatValues(handles);
% Update handles structure
guidata(hObject, handles);
% Bring main window to front
figure(handles.figure1)

% Fix stupid ginput
function [out1,out2,out3] = ginput(arg1)
out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
    if nargout == 1,
        if nargin == 1,
            out1 = trmginput(arg1);
        else
            out1 = trmginput;
        end
    elseif nargout == 2 || nargout == 0,
        if nargin == 1,
            [out1,out2] = trmginput(arg1);
        else
            [out1,out2] = trmginput;
        end
        if  nargout == 0
            out1 = [ out1 out2 ];
        end
    elseif nargout == 3,
        if nargin == 1,
            [out1,out2,out3] = trmginput(arg1);
        else
            [out1,out2,out3] = trmginput;
        end
    end
else
    
    fig = gcf;
    figure(gcf);
    
    if nargin == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  ischar(how_many) ...
                || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
                || ~(fix(how_many) == how_many) ...
                || how_many < 0
            error(message('MATLAB:ginput:NeedPositiveInt'))
        end
        if how_many == 0
            % If input argument is equal to zero points,
            % give a warning and return empty for the outputs.
            
            warning (message('MATLAB:ginput:InputArgumentZero'));
        end
    end
    
    % Setup the figure to disable interactive modes and activate pointers. 
    initialState = setupFcn(fig);
    
    % onCleanup object to restore everything to original state in event of
    % completion, closing of figure errors or ctrl+c. 
    c = onCleanup(@() restoreFcn(initialState));
       
    
    % We need to pump the event queue on unix
    % before calling WAITFORBUTTONPRESS
    drawnow
    char = 0;
    
    while how_many ~= 0
        % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
        try
            keydown = wfbp;
        catch %#ok<CTCH>
            waserr = 1;
        end
        if(waserr == 1)
            if(ishghandle(fig))
                cleanup(c);
                error(message('MATLAB:ginput:Interrupted'));
            else
                cleanup(c);
                error(message('MATLAB:ginput:FigureDeletionPause'));
            end
        end
        % g467403 - ginput failed to discern clicks/keypresses on the figure it was
        % registered to operate on and any other open figures whose handle
        % visibility were set to off
        figchildren = allchild(0);
        if ~isempty(figchildren)
            ptr_fig = figchildren(1);
        else
            error(message('MATLAB:ginput:FigureUnavailable'));
        end
        %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
        %         clicked figure has handlevisibility set to callback
        if(ptr_fig == fig)
            if keydown
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
            else
                button = get(fig, 'SelectionType');
                if strcmp(button,'open')
                    button = 1;
                elseif strcmp(button,'normal')
                    button = 1;
                elseif strcmp(button,'extend')
                    button = 2;
                elseif strcmp(button,'alt')
                    button = 3;
                else
                    error(message('MATLAB:ginput:InvalidSelection'))
                end
            end
            axes_handle = gca;
            drawnow;
            pt = get(axes_handle, 'CurrentPoint');
            
            how_many = how_many - 1;
            
            if(char == 13) % & how_many ~= 0)
                break;
            end
            
            out1 = [out1;pt(1,1)]; %#ok<AGROW>
            y = [y;pt(1,2)]; %#ok<AGROW>
            b = [b;button]; %#ok<AGROW>
        end
    end
    
    % Cleanup and Restore 
    cleanup(c);
    
    if nargout > 1
        out2 = y;
        if nargout > 2
            out3 = b;
        end
    else
        out1 = [out1 y];
    end
end    
    function key = wfbp
fig = gcf;
current_char = []; %#ok<NASGU>
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error(message('MATLAB:ginput:Interrupted'));
end
if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialState = setupFcn(fig)
% Store Figure Handle. 
initialState.figureHandle = fig; 
% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);
% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end
% Setup FullCrosshair Pointer without warning. 
oldwarnstate = warning('off', 'MATLAB:hg:Figure:Pointer');
% set(fig,'Pointer','fullcrosshair');
warning(oldwarnstate);
% Adding this to enable automatic updating of currentpoint on the figure 
set(fig,'WindowButtonMotionFcn',@(o,e) dummy());
% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    set(initialState.figureHandle,'WindowButtonMotionFcn','');
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
function dummy()
% do nothing, this is there to update the GINPUT WindowButtonMotionFcn. 
function cleanup(c)
if isvalid(c)
    delete(c);
end


function pushbuttonAbout_Callback(hObject, eventdata, handles)
msgbox(['These scripts were created by Colin Ophus [clophus@lbl.gov] ' ...
    'at the National Center for Electron Microscopy, ' ...
    'Lawrence Berkeley National Laboratory, ' ...
    'Berkeley, CA, USA.  Use at your own risk, ' ...
    'as accuracy cannot be guaranteed!']);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
