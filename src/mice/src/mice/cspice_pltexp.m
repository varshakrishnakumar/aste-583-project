function [overts] = cspice_pltexp( iverts, delta)
   switch nargin
      case 2
         iverts = zzmice_dp(iverts);
         delta  = zzmice_dp(delta);
      otherwise
         error ( [ 'Usage: [overts(3,3)] = ' ...
                   'cspice_pltexp( iverts(3,3), delta)' ] )
   end
   try
      [overts] = mice('pltexp_c', iverts, delta );
   catch spiceerr
      rethrow(spiceerr)
   end