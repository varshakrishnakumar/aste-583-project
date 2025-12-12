function [ npoint, alt ] = cspice_nearpt( positn, a, b, c )
   switch nargin
      case 4
         positn = zzmice_dp(positn);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
      otherwise
         error ( ['Usage: [_npoint(3)_, _alt_] = ' ...
                  'cspice_nearpt( _positn(3)_, a, b, c )'] )
   end
   try
      [nearpt] = mice( 'nearpt_s', positn, a, b, c  );
      npoint   = reshape( [nearpt.pos], 3, [] );
      alt      = reshape( [nearpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end