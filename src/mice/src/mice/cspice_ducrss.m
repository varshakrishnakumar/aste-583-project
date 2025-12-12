function [sout] = cspice_ducrss(s1, s2)
   switch nargin
      case 2
         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);
      otherwise
         error ( 'Usage: [_sout(6)_] = cspice_ducrss(_s1(6)_, _s2(6)_)' )
   end
   try
      [sout] = mice('ducrss_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end