function [xform] = cspice_twovxf( axdef, indexa, plndef, indexp )
   switch nargin
      case 4
         axdef  = zzmice_dp(axdef);
         indexa = zzmice_int(indexa);
         plndef = zzmice_dp(plndef);
         indexp = zzmice_int(indexp);
      otherwise
         error ( [ 'Usage: [xform(6,6)] = '                                 ...
                   'cspice_twovxf( axdef(6), indexa, plndef(6), indexp )' ] )
   end
   try
      [xform] = mice('twovxf_c', axdef, indexa, plndef, indexp);
   catch spiceerr
      rethrow(spiceerr)
   end