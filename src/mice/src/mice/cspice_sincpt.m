function [spoint, trgepc,                                               ...
          srfvec, found] = cspice_sincpt( method, target, et,   fixref, ...
                                          abcorr, obsrvr, dref, dvec  )
   switch nargin
      case 8
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);
      otherwise
         error( [ 'Usage: [ spoint, trgepc, srfvec, found] =  ' ...
                  'cspice_sincpt( `method`, `target`, et, `fixref`, ' ...
                                 '`abcorr`, `obsrvr`, `dref`, dvec)' ]  )
   end
   try
      [sincpt] = mice('sincpt_s', method, target, ...
                                  et, fixref, abcorr, obsrvr, dref, dvec);
      spoint = reshape( [sincpt.spoint], 3, [] );
      trgepc = reshape( [sincpt.trgepc], 1, [] );
      srfvec = reshape( [sincpt.srfvec], 3, [] );
      found  = reshape( [sincpt.found] , 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end