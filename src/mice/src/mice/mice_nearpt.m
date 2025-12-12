function [ npoint ] = mice_nearpt( positn, a, b, c )
   switch nargin
      case 4
         positn = zzmice_dp(positn);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
      otherwise
         error ( ['Usage: [_npoint_] = '                                   ...
                  'mice_nearpt( _positn(3)_, a, b, c )'] )
   end
   try
      [npoint] = mice( 'nearpt_s', positn, a, b, c  );
   catch spiceerr
      rethrow(spiceerr)
   end