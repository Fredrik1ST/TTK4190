% The function is copied from the free matlab package 'flypath3d' 
% http://www.wbint.pl/flypath3d/index.php
%
% Lines 93-97 were added in August 2023 (bogdan.l.hansen@ntnu.no)
% 
% NEW_OBJECT (new_object.m)
%
% Creates 3d object data set for visualization
% Original code was last updated: 2016-02-08
%
% SYNTAX:
%
%    new_object(filename,matrix,varargin)
%
% PARAMETERS:
%
%    filename    : output set file name (*.mat)
%    matrix      : array of kinematic vectors data provided in the NED (North-East-Down) frame, which is a common
%                  frame of reference used in aerospace and navigation applications. Ensure that the inputs are
%                  correctly transformed to NED coordinates before calling this function. 
%                  matrix(:,1) - x Cartesian coordinate vector [m],
%                  matrix(:,2) - y Cartesian coordinate vector [m],
%                  matrix(:,3) - z Cartesian coordinate vector [m],
%                  matrix(:,4) - roll angle vector [rad]
%                  matrix(:,5) - pitch angle vector [rad],
%                  matrix(:,6) - yaw angle vector [rad],
%
% OPTIONAL ARGUMENTS:
%
%    'alpha'     : alpha channel value (0-1 - default 1)
%    'edge'      : model edge color ([R G B] - default [.4 .4 .4])
%    'face'      : model face color ([R G B] - default [.5 .5 .5])
%    'model'     : 3d model file name (string - default 'missile.mat')
%    'path'      : trajectory line visibility on/off ('on','off' - default 'on')
%    'pathcolor' : trajectory line color ([R G B] - default [.3 .3 .3])
%    'pathwidth' : trajectory line width (default 1)
%    'scale'     : model scale (default 1)
%
% EXAMPLES OF USE:
%
%    trajectory = load('trajectory_tbm.mat');
%    new_object('tbm.mat',trajectory,...
%               'model','scud.mat','scale',5,...
%               'path','on','pathcolor',[.89 .0 .27]);
%

function new_object(filename,matrix,varargin)

  % --- check input parameters ---
  if nargin < 2
     error('Not enough input parameters!');
  end;
  
  % --- default input parameters ---
  pAlpha = 1.0;                 % alpha channel value: 0-1
  pEdge = [ .4 .4 .4 ];         % model edge color: [ R G B ]
  pFace = [ .5 .5 .5 ];         % model face color: [ R G B ]
  pModel = 'missile.mat';       % 3d model file name
  pPath = 'on';                 % trajectory line visibility
  pPathColor = [ .3 .3 .3 ];	% trajectory line color
  pPathWidth = 1;               % trajectory line width
  pScale = 1;                   % model scale

  % --- read input parameters ---
  i = 1;
  while i <= length(varargin)
     switch lower(varargin{i})
        case 'alpha'
           pAlpha = varargin{i+1};
		case 'edge'
		   pEdge = varargin{i+1};
		case 'face'
		   pFace = varargin{i+1};
		case 'model'
		   pModel = varargin{i+1};
		case 'path'
		   pPath = varargin{i+1};
		case 'pathcolor'
		   pPathColor = varargin{i+1};
        case 'pathwidth'
		   pPathWidth = varargin{i+1};
        case 'scale'
		   pScale = varargin{i+1};
         case 'enu'
           ENU = varargin{i+1};
     end
     i = i + 2;
  end;
  
  % --- kinematic value array ---
  pMatrix = matrix;
  
  % NED to ENU conversion necessary for visualization
  pMatrix = pMatrix(:,[2 1 3 4 5 6]);
  pMatrix(:,3) = -pMatrix(:,3);
  pMatrix(:,5) = -pMatrix(:,5);  
  pMatrix(:,6) = ssa(-pMatrix(:,6)+pi/2);
   
  % --- save data to file ---
  save(filename,'pMatrix',...
       'pAlpha','pEdge','pFace','pModel','pScale',...
       'pPath','pPathColor','pPathWidth');
end