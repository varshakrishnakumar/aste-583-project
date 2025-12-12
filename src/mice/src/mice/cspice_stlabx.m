function [corpos] = cspice_stlabx( pobj, vobs )
   switch nargin
      case 2
         pobj = zzmice_dp(pobj);
         vobs = zzmice_dp(vobs);
      otherwise
         error ( 'Usage: [corpos(3)] = cspice_stlabx( pobj(3), vobs(3) )' )
   end
   try
      [corpos] = mice('stlabx_c', pobj, vobs);
   catch spiceerr
      rethrow(spiceerr)
   end