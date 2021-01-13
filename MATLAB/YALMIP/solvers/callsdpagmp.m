function output = callsdpagmp(interfacedata)

% CALLSDPAGMP.m Call SDPA-GMP from YALMIP
% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    23/08/2016
%
%     Copyright (C) 2016  Giovanni Fantuzzi
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ----------------------------------------------------------------------- %

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Convert from internal (sedumi) format
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);

if options.verbose==0
    options.sdpa_gmp.print = 'no';
else
    options.sdpa_gmp.print = 'display';
end

if options.savedebug
    ops = options.sdpa_gmp;
    save sdpa_gmpdebug mDIM nBLOCK bLOCKsTRUCT c F ops
end

if options.showprogress
    showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL SDPA-GMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solvertime = tic;
[objVal,x,X,Y,INFO] = sdpagmp(mDIM,nBLOCK,bLOCKsTRUCT,c,F,[],[],[],options.sdpa_gmp);
solvertime = toc(solvertime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here onwards, like in YALMIP native callsdpa
% Create variables in YALMIP internal format
Primal = x;

Dual = [];
for i = 1:length(Y)
    Dual = [Dual;Y{i}(:)];
end

Slack = [];
if options.saveduals
    for i = 1:length(X)
        Slack = [Slack;X{i}(:)];
    end
end

switch (INFO.phasevalue)
    case 'pdOPT'
        problem = 0;
    case {'noINFO','pFEAS','dFEAS'}
        problem = 3;
    case {'pdFEAS'}
        problem = 4;
    case 'pFEAS_dINF'
        problem = 2;
    case 'pINF_dFEAS'
        problem = 1;
    case 'pUNBD'
        problem = 2;
    case 'dUNBD'
        problem = 1;
    case 'pdINF'
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
    solveroutput.objVal = objVal;
    solveroutput.x = x;
    solveroutput.X = X;
    solveroutput.Y = Y;
    solveroutput.INFO = INFO;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.mDIM = mDIM;
    solverinput.nBLOCK=nBLOCK;
    solverinput.bLOCKsTRUCT=bLOCKsTRUCT;
    solverinput.c=c;
    solverinput.F=F;
else
    solverinput = [];
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);


























