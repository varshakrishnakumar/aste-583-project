function [pnear, dist] = cspice_npedln( a, b, c, linept, linedr )
   switch nargin
      case 5
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
         linept = zzmice_dp(linept);
         linedr = zzmice_dp(linedr);
      otherwise
         error ( ['Usage: [ pnear(3), dist] = ' ...
                  'cspice_npedln( a, b, c, linept(3), linedr(3) )'] )
   end
   try
      [npedln] = mice( 'npedln_s',a, b, c, linept, linedr );
      pnear    = reshape( [npedln.pos], 3, [] );
      dist     = reshape( [npedln.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end