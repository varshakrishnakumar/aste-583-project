function [dvdot] = cspice_dvdot(s1, s2)
   switch nargin
      case 2
         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);
      otherwise
         error ( 'Usage: [_dvdot_] = cspice_dvdot(_s1(6)_, _s2(6)_)' )
   end
   try
      [dvdot] = mice('dvdot_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end