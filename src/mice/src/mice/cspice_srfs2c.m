function [code,found] = cspice_srfs2c( srfstr, bodstr)
   switch nargin
      case 2
         srfstr = zzmice_str(srfstr);
         bodstr = zzmice_str(bodstr);
      otherwise
         error ( 'Usage: [ code, found] = cspice_srfs2c( `srfstr`, `bodstr`)' )
   end
   try
      [srfs2c] = mice('srfs2c_s', srfstr, bodstr);
      code     = reshape( [srfs2c.code],  1, [] );
      found    = reshape( [srfs2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end