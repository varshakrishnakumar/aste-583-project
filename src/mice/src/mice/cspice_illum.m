function [phase, solar, emissn] = cspice_illum( target, et, abcorr, ...
                                                obsrvr, spoint )
   switch nargin
      case 5
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);
      otherwise
         error ( ['Usage: [_phase_, _solar_, _emissn_] = '           ...
                          'cspice_illum( `target`, _et_, `abcorr`, ' ...
                          '`obsrvr`, _spoint(3)_)'] )
   end
   try
      [phase, solar, emissn] = mice('illum_c', target, et, ...
                                               abcorr, obsrvr, spoint );
   catch spiceerr
      rethrow(spiceerr)
   end