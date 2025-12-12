function [rotate] = cspice_pxfrm2(from, to, etfrom, etto)
   switch nargin
      case 4
         from   = zzmice_str(from);
         to     = zzmice_str(to);
         etfrom = zzmice_dp(etfrom);
         etto   = zzmice_dp(etto);
      otherwise
         error ( [ 'Usage: [_rotate(3,3)_] = ' ...
                   'cspice_pxfrm2( `from`, `to`, _etfrom_, _etto_ )' ] )
   end
   try
      [rotate] = mice('pxfrm2_c',from,to,etfrom,etto);
   catch spiceerr
      rethrow(spiceerr)
   end