function [ limb ] = cspice_edlimb( a, b, c, viewpt )
   switch nargin
      case 4
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
         viewpt = zzmice_dp(viewpt);
      otherwise
         error ( 'Usage: [ limb ] = cspice_edlimb( a, b, c, viewpt )' )
   end
   try
      [ limb ] = mice( 'edlimb_s', a, b, c, viewpt );
   catch spiceerr
      rethrow(spiceerr)
   end