function [state, lt] = cspice_spkcvt( trgsta, trgepc, trgctr, ...
                                      trgref, et,     outref, ...
                                      refloc, abcorr, obsrvr )
   switch nargin
      case 9
         trgsta = zzmice_dp(trgsta);
         trgepc = zzmice_dp(trgepc);
         trgctr = zzmice_str(trgctr);
         trgref = zzmice_str(trgref);
         et     = zzmice_dp(et);
         outref = zzmice_str(outref);
         refloc = zzmice_str(refloc);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error( ['Usage: [state(6), lt] = cspice_spkcvt( trgsta(6), ' ...
                         'trgepc, `trgctr`, `trgref`, et, `outref`, ' ...
                         '`refloc`, `abcorr`, `obsrvr` )'] )
   end
   try
      [starg] = mice('spkcvt_s', trgsta, trgepc, trgctr, ...
                                 trgref, et,     outref, ...
                                 refloc, abcorr, obsrvr );
      state   = reshape( [starg.state], 6, [] );
      lt      = reshape( [starg.lt   ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end