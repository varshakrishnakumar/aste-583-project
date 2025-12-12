function [unitim] = cspice_unitim( epoch, insys, outsys )
   switch nargin
      case 3
         epoch  = zzmice_dp(epoch);
         insys  = zzmice_str(insys);
         outsys = zzmice_str(outsys);
      otherwise
         error( [ 'Usage: [_unitim_] = '                                    ...
                  'cspice_unitim( _epoch_, `insys`, `outsys` )' ] )
   end
   try
      [unitim] = mice( 'unitim_c', epoch, insys, outsys);
   catch spiceerr
      rethrow(spiceerr)
   end