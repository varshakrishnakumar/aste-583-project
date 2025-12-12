function [vdist] = cspice_vdist( v1, v2 )
   switch nargin
      case 2
         v1 = zzmice_dp(v1);
         v2 = zzmice_dp(v2);
      otherwise
         error ( 'Usage: [_vdist_] = cspice_vdist( _v1(3)_, _v2(3)_ )' )
   end
   try
      [vdist] = mice('vdist_c',v1, v2);
   catch spiceerr
      rethrow(spiceerr)
   end