function [handle, descr, ident, found] = cspice_spksfs(body, et)
   switch nargin
      case 2
         body = zzmice_int(body);
         et   = zzmice_dp(et);
      otherwise
         error ( ['Usage: [handle, descr(5), `ident`, found] =' ...
                              ' cspice_spksfs( body, et)'] )
   end
   try
      spksfs = mice( 'spksfs_s', body, et );
      handle  = reshape( [spksfs.handle], 1, [] );
      descr   = reshape( [spksfs.descr ], 5, [] );
      ident   = char( spksfs.ident );
      found   = reshape( [spksfs.found ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end