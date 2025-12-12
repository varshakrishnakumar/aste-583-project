function [plane] = cspice_nvp2pl( normal, point )
   switch nargin
      case 2
         normal = zzmice_dp(normal);
         point  = zzmice_dp(point);
      otherwise
         error ( ['Usage: [plane] = cspice_nvp2pl( normal(3), point(3) )'] )
   end
   try
      [plane] = mice('nvp2pl_c', normal, point );
   catch spiceerr
      rethrow(spiceerr)
   end