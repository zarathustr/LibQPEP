function install_sdpt3( varargin )

%SDPT3 installation script
%
% Heavy rewrite by Michael C. Grant, pulled from similar work on SeDuMi.
%

need_rebuild = any( strcmp( varargin, '-rebuild' ) );
no_path = any( strcmp( varargin, '-nopath' ) );

targets64={...
    'mexMatvec', 'mexProd2', 'mexProd2nz', 'mexexpand', 'mexinprod', ...
    'mexmat', 'mexnnz', 'mexqops', 'mexschur', 'mexschurfun', ...
    'mexskron', 'mexsmat', 'mexsvec', 'mextriang', 'mextriangsp' ...
    };

fs = filesep;
mpath = mfilename('fullpath');
mpath = mpath( 1 : max([1,strfind(mpath,fs)]) - 1 );
mbase = [ mpath, fs, 'Solver', fs, 'Mexfun' ];
mdir = '';
ISOCTAVE = exist('OCTAVE_VERSION','builtin');
VERSION = [1,0.01]*sscanf(version,'%d.%d');
if ISOCTAVE, prog = 'Octave'; else prog = 'Matlab'; end
COMPUTER = computer;
mext = mexext;

%
% We don't want to rebuild the binaries if they're already present, unless
% the user has specifically asked for a rebuild. Frankly, we can't
% guarantee that rebuilding will work.
%

if ISOCTAVE,
    page_output_immediately( true, 'local' );
end

line = '---------------------------------------------------------------------------';
fprintf( '\n%s\nSDPT3 installation script\n   Directory: %s\n   %s %s on %s\n%s\n', ...
    line, mpath, prog, version, COMPUTER, line );

if ~need_rebuild,
    fprintf( 'Looking for existing binaries...' );
    mdir = '';
    if ISOCTAVE && VERSION >= 4,
        switch computer,
        case 'i686-w64-mingw32',
            mdir = 'o_win';
        end
    end
    if ~isempty(mdir) && ~exist( [ mbase, fs, mdir ], 'dir' ),
        mdir = '';
    end
    nfound = [ 0, 0 ];
    for k = 1 : length(targets64),
        targ = [ targets64{k}, '.', mext ];
        if exist( [ mbase, fs, targ ], 'file' ),
            nfound(1) = nfound(1) + 1;
        elseif ~isempty(mdir) && exist( [ mbase, fs, mdir, fs, targ ], 'file' ),
            nfound(2) = nfound(2) + 1;
        end
    end
    if sum(nfound) == 0,
        fprintf( 'none found; building...\n' );
        need_rebuild = true;
    elseif sum(nfound) < length(targets64),
        fprintf( 'incomplete set found.\n' );
        disp( line );
        warning(['Some of the binaries for this platform were found, but some', char(10), ...
                 'were missing as well. This may mean your download was corrupt;', char(10), ...
                 'consider downloading and unpacking SDPT3 again. Otherwise, to', char(10), ...
                 'try rebuilding the MEX files yourself, run this command:', char(10), ...
                 '    install_sdpt3 -rebuild', char(10), ...
                 line, char(10), char(10)]);
        return;
    else
        fprintf( 'found!\n' );
        fprintf( '   If for some reason you need to rebuild the binaries, use this command:\n' );
        fprintf( '      install_sdpt3 -rebuild\n' );
    end
else
    nfound = [1,0];
end

if need_rebuild,
    disp( 'Attempting to recompile the SDPT3 binaries:' );
    % Note the use of 0.01 here. That's because version 7 had more than 10
    % minor releases, so 7.10-7.14 need to be ordered after 7.01-7.09.
    libs = {};
    flags = {'-O'};
    if ISOCTAVE,
        % Octave has mwSize and mwIndex hardcoded in mex.h as ints.
        % There is no definition for mwSignedIndex so include it here.  
        % This means that Octave ignores the -largeArrayDims flag.
        if VERSION < 3.08,
            flags{end+1} = '-DmwSignedIndex=int';
        end
    else
        if ispc,
            flags = {'-DPC'};
        elseif isunix,
            flags = {'-DUNIX'};
        end
        if strcmp(COMPUTER(end-1:end),'64') && ( VERSION >= 7.03 ),
            flags{end+1} = '-largeArrayDims';
        elseif VERSION < 7.03,
            flags{end+1} = '-DmwIndex=int';
            flags{end+1} = '-DmwSize=int';
            flags{end+1} = '-DmwSignedIndex=int';
        end
        if VERSION >= 7 && ispc,
            if IS64BIT, dirval = 'win64'; else dirval = 'win32'; end
            libdir = [ matlabroot, fs, 'extern', fs, 'lib', fs, dirval, fs ];
            if exist( [ libdir, 'microsoft' ], 'dir' ),
                libdir = [ libdir, 'microsoft' ];
                found = true;
            elseif exist( [ libdir, 'msvc60' ], 'dir' ),
                libdir = [ libdir, 'msvc60' ];
                found = true;
            elseif exist( [ libdir, 'lcc' ], 'dir' ),
                libdir = [ libdir, 'lcc' ];
                found = true;
            end
            if found,
                libs{end+1} = [ '-L"', libdir, '"' ];
            end
        end
        if VERSION >= 7.05,
            libs{end+1} = '-lmwblas';
        else
            libs{end+1} = '-lmwlapack';
        end
    end
    libs = sprintf( ' %s', libs{:} );
    flags = sprintf( ' %s', flags{:} );
    olddir = pwd;
    cd( mbase );
    failed = false;
    fprintf( 'Template: mex%s <sources>%s\n', flags, libs );
    for i=1:length(targets64),
        targ = targets64{i};
        mfile = [ targ(1:min(strfind(targ,'.'))), mext ];
        temp = [ 'mex ', flags, ' ', targets64{i}, '.c ', libs ];
        fprintf( '   %s: %s\n', mfile, targ );
        eval( temp, 'failed=true;' ); %#ok
    end
    cd( olddir );
    if failed,
        fprintf( 'At least one compilation failure occurred.\n' ); %#ok
        nfound = [0,0];
    else
        fprintf( 'Compilation successful.\n' );
        nfound = [1,0];
    end
end

if ~any(nfound),
    disp( line );
    error(['SDPT3 was not successfully installed.', char(10), ...
    	   'Please attempt to correct the errors and try again.']);
elseif ~no_path,
    disp( line );
    fprintf( 'Adding SDPT3 to the %s path:\n', prog );
    paths = { 'Base', 'Solver', 'HSDSolver', 'Binaries', 'Binaries', 'Examples' ; ...
              '', 'Solver', 'HSDSolver', [ 'Solver', fs, 'Mexfun' ], ...
              [ 'Solver', fs, 'Mexfun', fs, mdir ], 'Examples' ; ...
              0, 0, 0, 1, 2, 0 };
    ps = pathsep;
    pp = [ ps, path, ps ];
    already = true;
    for k = 1 : size(paths,2),
        if paths{3,k} ~= 0 && nfound(paths{3,k}) == 0, continue; end
        fprintf( '   %s...', paths{1,k} );
        t_dir = mpath;
        if ~isempty(paths{2,k}), 
            t_dir = [ t_dir, fs, paths{2,k} ]; %#ok
        end
        if ~any(strfind(pp,[ps,t_dir,ps])),
            already = false;
            pp = [ pp, t_dir, ps ]; %#ok
            fprintf( 'added.\n' );
        else
            fprintf( 'already there.\n' );
        end
    end
    if ~already,
        path(pp);
        fprintf( 'Please save the %s path if you want to use SDPT3 from any directory.\n', prog );
    end
    disp( line );
    disp('SDPT3 has been succesfully installed.' );
    disp( 'For more information, type "help sdpt3" or see the user guide.')
end

fprintf('%s\n\n',line);
