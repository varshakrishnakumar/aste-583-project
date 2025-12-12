function tictoc()
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   cont        = SPICETRUE;
   prec        = 3;
   utc = strvcat( '9 JAN 1986 03:12:59.22451', ...
                  '1/9/86 3:12:59.22451',      ...
                  '86-365//12:00',             ...
                  'JD 2451545',                ...
                  '77 JUL 1',                  ...
                  '1 JUL ''29'                 ...
                  );
   disp( ' '                                                         )
   disp( '                 Welcome to TICTOC'                        )
   disp( ' '                                                         )
   disp( 'This program demonstrates the use of the time conversion ' )
   disp( 'utility routines: cspice_str2et and cspice_et2utc.'        )
   disp( ' '                                                         )
   disp( ' ' )
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )
   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )
   et = cspice_str2et( utc );
   NCASES = size( et, 2 );
   i      = 1;
   while ( (i <= NCASES)  &&  (cont == SPICETRUE) )
      txt = sprintf( '      Example UTC time      :  %s', utc(i,:) );
      disp( txt)
      disp( ' ' )
      txt = sprintf( '      Corresponding ET      :  %f', et(i)    );
      disp( txt)
      format  = 'C';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC calendar format   :  %s', timestr );
      disp( txt)
      format  = 'D';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC day of year format:  %s', timestr );
      disp( txt)
      format  = 'J';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC day of year format:  %s', timestr );
      disp( txt)
      if ( i < NCASES )
         disp( ' ' )
         answer = input( 'Continue? (Enter Y or N): ', 's');
         disp( ' ' )
      end
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end
      i  = i + 1;
   end
   cspice_kclear