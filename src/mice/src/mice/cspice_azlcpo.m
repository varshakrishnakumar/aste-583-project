function [azlsta, lt] = cspice_azlcpo( method, target, et,                 ...
                                       abcorr, azccw,  elplsz,             ...
                                       obspos, obsctr, obsref )
   switch nargin
      case 9
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);
         obspos = zzmice_dp(obspos);
         obsctr = zzmice_str(obsctr);
         obsref = zzmice_str(obsref);
      otherwise
         error ( [ 'Usage: [azlsta(6), lt] = '                              ...
                   'cspice_azlcpo( `method`, `target`, et, `abcorr`, '      ...
                   'azccw, elplsz, obspos(3), `obsctr`, `obsref` )' ] )
   end
   try
      [azlsta, lt] = mice('azlcpo_c', method, target, et,     abcorr,       ...
                          azccw,      elplsz, obspos, obsctr, obsref);
   catch spiceerr
      rethrow(spiceerr)
   end