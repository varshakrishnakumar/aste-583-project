function [ellips] = cspice_cgv2el( center, vec1, vec2 )
   switch nargin
      case 3
         center = zzmice_dp(center);
         vec1   = zzmice_dp(vec1);
         vec2   = zzmice_dp(vec2);
      otherwise
         error ( ['Usage: [ellips] = ' ...
                  'cspice_cgv2el( center(3), vec1(3), vec2(3) )'] )
   end
   try
      [ellips] = mice('cgv2el_c', center, vec1, vec2 );
   catch spiceerr
      rethrow(spiceerr)
   end