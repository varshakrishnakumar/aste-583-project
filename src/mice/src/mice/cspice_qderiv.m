function [dfdt] = cspice_qderiv( f0, f2, delta )
   switch nargin
      case 3
         f0    = zzmice_dp(f0);
         f2    = zzmice_dp(f2);
         delta = zzmice_dp(delta);
      otherwise
         error ( 'Usage: [dfdt(N)] = cspice_qderiv( f0(N), f2(N), delta )' )
   end
   try
      [dfdt] = mice('qderiv_c', f0, f2, delta);
   catch spiceerr
      rethrow(spiceerr)
   end