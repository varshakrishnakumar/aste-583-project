function [plane] = cspice_psv2pl( point, span1, span2 )
   switch nargin
      case 3
         point = zzmice_dp(point);
         span1 = zzmice_dp(span1);
         span2 = zzmice_dp(span2);
      otherwise
         error ( ['Usage: [plane] = ' ...
                  'cspice_psv2pl( point(3), span1(3), span2(3) )'] )
   end
   try
      [plane] = mice('psv2pl_c', point, span1, span2 );
   catch spiceerr
      rethrow(spiceerr)
   end