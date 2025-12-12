function [state, lt] = cspice_spkcpt( trgpos, trgctr, trgref, ...
                                      et,     outref, refloc, ...
                                      abcorr, obsrvr )
   switch nargin
      case 8
         trgpos = zzmice_dp(trgpos);
         trgctr = zzmice_str(trgctr);
         trgref = zzmice_str(trgref);
         et     = zzmice_dp(et);
         outref = zzmice_str(outref);
         refloc = zzmice_str(refloc);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error( ['Usage: [ state(6), lt] = cspice_spkcpt( ' ...
                          'trgpos(3), `trgctr`, `trgref`, ' ...
                          'et, `outref`, `refloc`, '        ...
                          '`abcorr`, `obsrvr` )'] )
   end
   try
      [starg] = mice('spkcpt_s', trgpos, trgctr, trgref, ...
                                 et,     outref, refloc, ...
                                 abcorr, obsrvr );
      state   = reshape( [starg.state], 6, [] );
      lt      = reshape( [starg.lt   ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end