function [exedir,dlgshown]=mbe_getprgmdir(prefstr)

% e.g., prefstr='sweepfinderrun_prgmdir';
exedir='';
dlgshown=false;
dirok=false;
if ispref('mbetoolbox',prefstr)
    exedir=getpref('mbetoolbox',prefstr);
    if isdir(exedir), dirok=true; end
end
if ~dirok
        dlgshown=true;
    	ButtonName=questdlg(sprintf('Do you want to set up external program path?'),...
        			    'External Program First-time Run',...
			            'Yes','No','Cancel','Yes');
        if strcmp(ButtonName,'Yes')
           exedir=uigetdir(pwd);
                    if exedir~=0
                       setpref('mbetoolbox',prefstr,exedir);
                       helpdlg('You are ready to run.','External Program Path Saved')
                    end
        end
        return;
end
