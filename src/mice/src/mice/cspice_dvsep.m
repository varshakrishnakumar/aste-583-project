function [dvsep] = cspice_dvsep(s1, s2)
   switch nargin
      case 2
         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);
      otherwise
         error ( 'Usage: [_dvsep_] = cspice_dvsep(_s1(6)_, _s2(6)_)' )
   end
   try
      [dvsep] = mice('dvsep_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end