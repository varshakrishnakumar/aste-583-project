function [calstr] = cspice_etcal(et)
   switch nargin
      case 1
         et = zzmice_dp(et);
      otherwise
         error( 'Usage: [_`calstr`_] = cspice_etcal(_et_)' )
   end
   try
      [calstr] = mice('etcal_c', et );
   catch spiceerr
      rethrow(spiceerr)
   end