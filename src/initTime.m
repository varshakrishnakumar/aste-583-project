function [ET0, LTM, LCM] = initTime()
    persistent spiceLoaded
    if isempty(spiceLoaded) || ~spiceLoaded
        thisDir = fileparts(mfilename('fullpath'));
        miceDir = fullfile(thisDir, 'mice');
        addpath(genpath(fullfile(miceDir,'src','mice')));
        addpath(genpath(fullfile(miceDir,'lib')));
        addpath(genpath(fullfile(miceDir,'include')));
        addpath(genpath(fullfile(miceDir,'generic_kernels')));
        cspice_kclear;
        kernDir = fullfile(miceDir,'generic_kernels');
        lsk = "";
        lskCandidates = ["naif0012.tls", "naif0012.tls.pc", "naif0012.tls.pc.txt"];
        for k = lskCandidates
            f = fullfile(kernDir, k);
            if exist(f, 'file')
                lsk = f; break;
            end
        end
        if lsk == ""
            error('initTime:MissingLSK', 'No leapseconds kernel found in %s', kernDir);
        end
        cspice_furnsh(char(lsk));
        de = fullfile(kernDir, 'de440s.bsp');
        if ~exist(de,'file')
            error('initTime:MissingDE', 'Missing ephemeris kernel %s', de);
        end
        cspice_furnsh(de);
        spiceLoaded = true;
    end
    ET0 = cspice_str2et('2025 DEC 01 00:00:00 UTC');
    LTM = cspice_str2et('2025 DEC 16 00:00:00 UTC');
    LCM = cspice_str2et('2025 DEC 25 00:00:00 UTC');
end