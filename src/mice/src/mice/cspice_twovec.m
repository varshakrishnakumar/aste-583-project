function [mout] = cspice_twovec(axdef, indexa, plndef, indexp)
   switch nargin
      case 4
         axdef  = zzmice_dp(axdef);
         indexa = zzmice_int(indexa);
         indexp = zzmice_int(indexp);
         plndef = zzmice_dp(plndef);
      otherwise
         error ( ['Usage: [mout(3,3)] = cspice_twovec( axdef(3), ' ...
                                       'indexa, plndef(3), indexp)'] )
   end
   try
      [mout] = mice('twovec_c', axdef, indexa, plndef, indexp);
   catch spiceerr
      rethrow(spiceerr)
   end