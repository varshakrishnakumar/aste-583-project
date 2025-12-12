function [visibl] = cspice_fovtrg( inst,   target, tshape, ...
                                   tframe, abcorr, obsrvr, et )
    switch nargin
        case 7
            inst   = zzmice_str(inst);
            target = zzmice_str(target);
            tshape = zzmice_str(tshape);
            tframe = zzmice_str(tframe);
            abcorr = zzmice_str(abcorr);
            obsrvr = zzmice_str(obsrvr);
            et     = zzmice_dp(et);
        otherwise
            error ( ['Usage: [_visibl_] = ' ...
                  'cspice_fovtrg( `inst`, `target`, '    ...
                  '`tshape`, `tframe`, `abcorr`, ' ...
                  '`obsrvr`, _et_]' ] )
   end
   try
      [visibl] = mice('fovtrg_c', inst,   target, tshape, ...
                                  tframe, abcorr, obsrvr, et );
      [visibl] = zzmice_logical( visibl );
   catch spiceerr
      rethrow(spiceerr)
   end