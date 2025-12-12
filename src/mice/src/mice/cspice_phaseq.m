function [phaseq] = cspice_phaseq( et, target, illmn, obsrvr, abcorr )
   switch nargin
      case 5
         et     = zzmice_dp(et);
         target = zzmice_str(target);
         illmn = zzmice_str(illmn);
         obsrvr = zzmice_str(obsrvr);
         abcorr = zzmice_str(abcorr);
      otherwise
         error ( ['Usage: [_phaseq_] = cspice_phaseq( _et_, '       ...
                                  '`target`, `illmn`, `obsrvr`, `abcorr` )'] )
   end
   try
      [phaseq] = mice('phaseq_c', et, target, illmn, obsrvr, abcorr );
   catch spiceerr
      rethrow(spiceerr)
   end