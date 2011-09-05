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

function run
%function run

global fig2;

dir = cd; %gets current directory
addpath([dir '\model']); %looks for functions in the sub-folder .\model
addpath([dir '\GUI']); %looks for functions in the sub-folder .\GUI

fig2 = GUI_ControlPanel2D;
handles = guihandles(fig2);
GUI_ControlPanel2D('resetbutton_Callback',fig2,0,handles);

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010,2011. Guillaume Riflet,
%Instituto Superior Técnico da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------