function [sout] = cspice_dvhat( s1 )
   switch nargin
      case 1
         s1 = zzmice_dp(s1);
      otherwise
         error ( 'Usage: [_sout(6)_] = cspice_dvhat( _s1(6)_ )' )
   end
   try
      [sout] = mice('dvhat_c',s1);
   catch spiceerr
      rethrow(spiceerr)
   end