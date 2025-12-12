function [spoint, dist, trgepc, obspos, found] = ...
          cspice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)
   switch nargin
      case 7
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);
      otherwise
         error( [ 'Usage: [_spoint(3)_, _dist_, _trgepc_, '   ...
                  '_obspos(3)_, _found_ ] = '                 ...
                  'cspice_srfxpt( `method`, `target`, _et_, ' ...
                  '`abcorr`, `obsrvr`, `dref`, dvec(3))']  )
   end
   try
      [surf] = mice('srfxpt_s', method, target, et,   ...
                     abcorr,    obsrvr, dref,   dvec);
      spoint = reshape( [surf.spoint], 3, [] );
      dist   = reshape( [surf.dist]  , 1, [] );
      trgepc = reshape( [surf.trgepc], 1, [] );
      obspos = reshape( [surf.obspos], 3, [] );
      found  = reshape( [surf.found] , 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end