function [appobj] = cspice_stelab( pobj, vobs )
   switch nargin
      case 2
         pobj = zzmice_dp(pobj);
         vobs = zzmice_dp(vobs);
      otherwise
         error ( 'Usage: [appobj(3)] = cspice_stelab( pobj(3), vobs(3) )' )
   end
   try
      [appobj] = mice('stelab_c', pobj, vobs);
   catch spiceerr
      rethrow(spiceerr)
   end