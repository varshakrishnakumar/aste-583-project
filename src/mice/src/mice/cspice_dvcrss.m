function [sout] = cspice_dvcrss(s1, s2)
   switch nargin
      case 2
         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);
      otherwise
         error ( 'Usage: [_sout(6)_] = cspice_dvcrss(_s1(6)_, _s2(6)_)' )
   end
   try
      [sout] = mice( 'dvcrss_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end