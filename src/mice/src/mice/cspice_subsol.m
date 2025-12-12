function [spoint] = cspice_subsol( method, target, et, abcorr, obsrvr )
   switch nargin
      case 5
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_spoint(3)_ ] = '            ...
                  'cspice_subsol( `method`, `target`, ' ...
                                  '_et_, `abcorr`, `obsrvr`)']  )
   end
   try
      [spoint] = mice('subsol_c', method, target, et, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end