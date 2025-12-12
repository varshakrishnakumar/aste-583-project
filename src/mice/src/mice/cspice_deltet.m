function [delta] = cspice_deltet( epoch, eptype )
   switch nargin
      case 2
         eptype = zzmice_str(eptype);
         epoch  = zzmice_dp(epoch);
      otherwise
         error ( 'Usage: [_delta_] = cspice_deltet( _epoch_, `eptype`)' )
   end
   try
      [delta] = mice('deltet_c', epoch, eptype) ;
   catch spiceerr
      rethrow(spiceerr)
   end