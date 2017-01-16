function [output] = list2celldlg

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

   prompt={'List:'};
   def={''};
   dlgTitle='List2Cell';
   lineNo=10;
   AddOpts.Resize='off';
   AddOpts.WindowStyle='modal';
   AddOpts.Interpreter='none';
   answer=inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);
   
 if ~isempty(answer)
     output=cellstr(answer{1});
 else 
     output=[];
 end
