%
%	startup file
%

disp([' Running ' ...
      which('startup')]);

clear path;
dirstr = '';
% GENERIC STUFF (change these to your directory, unless you are on CS-lab)
global matlabVisRoot;
matlabVisRoot = [userpath_fix,filesep,'csc2503_hw',filesep,'matlabVisTools'];

%dirstr=[matlabVisRoot];  % Put the Matlab root on the path
rootstr=[pathsep, matlabVisRoot];  % prefix for subsequent directories

% IMAGES directory (sample images)
%dirstr=[dirstr,rootstr,'/images'];
dirstr = [matlabVisRoot,filesep,'images'];

% ISE TOOLBOX DIRECTORY, TEACHING VERSION
rootstr=[pathsep,matlabVisRoot,filesep,'iseToolbox'];
dirstr=[dirstr,rootstr,filesep,'color'];
dirstr=[dirstr,rootstr,filesep,'filters'];
dirstr=[dirstr,rootstr,filesep,'graphics'];
dirstr=[dirstr,rootstr,filesep,'imops'];
dirstr=[dirstr,rootstr,filesep,'matrix'];
dirstr=[dirstr,rootstr,filesep,'pyrTools',filesep,'MEX'];
dirstr=[dirstr,rootstr,filesep,'pyrTools'];
dirstr=[dirstr,rootstr,filesep,'stats'];
dirstr=[dirstr,rootstr,filesep,'synth'];
dirstr=[dirstr,rootstr,filesep,'tutorials'];
dirstr=[dirstr,rootstr,filesep,'utility'];

% UTVIS TOOLBOX DIRECTORY
rootstr=[pathsep,matlabVisRoot,filesep,'utvisToolbox'];
dirstr=[dirstr,rootstr,filesep,'colour'];
dirstr=[dirstr,rootstr,filesep,'edge'];
dirstr=[dirstr,rootstr,filesep,'eigenEyes'];
dirstr=[dirstr,rootstr,filesep,'file'];
dirstr=[dirstr,rootstr,filesep,'tutorials'];

path(pathdef);
path([dirstr,pathsep,path]);
clear p dirstr rootstr;
clear startup.m;
