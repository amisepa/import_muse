% eegplugin_import_muse() - EEGLAB plugin for importing Muse data recorded
%             with the Muse Monitor or Muse Direct Apps.
%
% Usage:
%   >> eegplugin_import_muse(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Cedric Cannard, CerCo, CNRS
%
% Copyright (C) 2021 Cedric Cannard
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function vers = eegplugin_import_muse(fig, trystrs, catchstrs)

vers = 'import_muse1.0';
if nargin < 3
    error('eegplugin_import_muse requires 3 arguments');
end

%Add folder to path
p = which('eegplugin_import_muse.m');
p = p(1:strfind(p,'eegplugin_import_muse.m')-1);
if ~exist('eegplugin_import_muse','dir')
    addpath(p);
end

%Find import data menu
menui = findobj(fig, 'tag', 'import data');

%Menu callbacks
% try
% comcnt = [trystrs.no_check '[EEGTMP LASTCOM] = import_muse;'  catchstrs.new_non_empty];
% catch
comcnt = [trystrs.no_check '[EEGTMP ACCTMP GYRTMP PPGTMP AUXTMP LASTCOM] = import_muse;'  catchstrs.new_non_empty];
% end

%Create menus
uimenu(menui, 'label', 'MUSE .csv file (from Mind Monitor or Muse Direct)', 'separator', 'on', 'callback', comcnt);

end
