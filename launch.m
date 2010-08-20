function launch
%function launch

global fig2;

fig2 = ControlPanel2D;
handles = guihandles(fig2);
ControlPanel2D('resetbutton_Callback',fig2,0,handles);