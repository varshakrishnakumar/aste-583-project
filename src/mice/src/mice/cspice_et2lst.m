function [hr, mn, sc, time, ampm] = cspice_et2lst( et, body, lon, type)
   switch nargin
      case 4
         et    = zzmice_dp(et);
         body  = zzmice_int(body);
         lon   = zzmice_dp(lon);
         type  = zzmice_str(type);
      otherwise
         error ( ['Usage: [ _hr_, _mn_, _sc_, _`time`_, _`ampm`_] = ' ...
                 'cspice_et2lst( _et_, body, lon, `type`)'] )
   end
   try
      [hr, mn, sc, time, ampm] = mice('et2lst_c', et, body, lon, type);
      hr = zzmice_dp(hr);
      mn = zzmice_dp(mn);
      sc = zzmice_dp(sc);
   catch spiceerr
      rethrow(spiceerr)
   end