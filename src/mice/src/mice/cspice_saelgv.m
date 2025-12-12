function [ smajor, sminor ] = cspice_saelgv( vec1, vec2 )
   switch nargin
      case 2
         vec1   = zzmice_dp(vec1);
         vec2   = zzmice_dp(vec2);
      otherwise
         error ( ['Usage: [ smajor(3), sminor(3) ] = ' ...
                  'cspice_saelgv( vec1(3), vec2(3) )'] )
   end
   try
      [ smajor, sminor ] = mice('saelgv_c', vec1, vec2 );
   catch spiceerr
      rethrow(spiceerr)
   end