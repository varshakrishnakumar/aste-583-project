function [visibl] = cspice_fovray( inst,   raydir, rframe,  ...
                                   abcorr, obsrvr, et )
    switch nargin
        case 6
            inst    = zzmice_str(inst);
            raydir  = zzmice_dp (raydir);
            rframe  = zzmice_str(rframe);
            abcorr  = zzmice_str(abcorr);
            obsrvr  = zzmice_str(obsrvr);
            et      = zzmice_dp(et);
        otherwise
            error ( ['Usage: [_visibl_] = '              ...
                  'cspice_xfmsta( `inst`, _raydir[6]_, ' ...
                  '`rframe`, `abcorr`, `obsrvr`, _et_]' ] )
   end
   try
      [visibl] = mice('fovray_c', inst, raydir,   rframe, ...
                                  abcorr,     obsrvr, et );
      [visibl] = zzmice_logical( visibl );
   catch spiceerr
      rethrow(spiceerr)
   end