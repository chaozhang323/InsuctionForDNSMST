	  fG  ¨   k820309    ë          12.0        %5eY                                                                                                           
       io.f90 DECOMP_2D_IO              gen@DECOMP_2D_WRITE_ONE gen@DECOMP_2D_READ_ONE gen@DECOMP_2D_WRITE_VAR gen@DECOMP_2D_READ_VAR gen@DECOMP_2D_WRITE_SCALAR gen@DECOMP_2D_READ_SCALAR gen@DECOMP_2D_WRITE_PLANE gen@DECOMP_2D_WRITE_EVERY                      @                             
                                                          
                                                              u #WRITE_ONE_REAL    #WRITE_ONE_COMPLEX 
   #         @     @X                                              #WRITE_ONE_REAL%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                 PRESENT           
   @                                                   
@ @@                                                
              &                   &                   &                                                     
@ @@                                                 1           
 @@                                   è             #DECOMP_INFO 	   #         @     @X                            
                  #WRITE_ONE_COMPLEX%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                 PRESENT           
   @                                                   
@ @@                                                              &                   &                   &                                                     
@ @@                                                 1           
 @@                                   è             #DECOMP_INFO 	                                                          u #READ_ONE_REAL    #READ_ONE_COMPLEX    #         @     @X                                              #READ_ONE_REAL%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                 PRESENT           
   @                                                   
D @@                                                
 	              &                   &                   &                                                     
@ @@                                                 1           
 @@                                   è             #DECOMP_INFO 	   #         @     @X                                              #READ_ONE_COMPLEX%PRESENT    #IPENCIL    #VAR    #FILENAME    #OPT_DECOMP                  @                                 PRESENT           
   @                                                   
D @@                                                               &                   &                   &                                                     
@ @@                                                 1           
 @@                                   è             #DECOMP_INFO 	                                                          u #WRITE_VAR_REAL    #WRITE_VAR_COMPLEX #   #         @     @X                                              #WRITE_VAR_REAL%PRESENT    #FH    #DISP    #IPENCIL     #VAR !   #OPT_DECOMP "                 @                                 PRESENT           
@ @@                                                   
D @@                                                   
   @                                                    
@ @@                             !                   
              &                   &                   &                                                     
 @@                              "     è             #DECOMP_INFO 	   #         @     @X                            #                  #WRITE_VAR_COMPLEX%PRESENT $   #FH %   #DISP &   #IPENCIL '   #VAR (   #OPT_DECOMP )                 @                            $     PRESENT           
@ @@                              %                     
D @@                             &                      
   @                              '                     
@ @@                             (                                 &                   &                   &                                                     
 @@                              )     è             #DECOMP_INFO 	                                                          u #READ_VAR_REAL *   #READ_VAR_COMPLEX 1   #         @     @X                            *                  #READ_VAR_REAL%PRESENT +   #FH ,   #DISP -   #IPENCIL .   #VAR /   #OPT_DECOMP 0                 @                            +     PRESENT           
@ @@                              ,                     
D @@                             -                      
   @                              .                     
D @@                             /                   
               &                   &                   &                                                     
 @@                              0     è             #DECOMP_INFO 	   #         @     @X                            1                  #READ_VAR_COMPLEX%PRESENT 2   #FH 3   #DISP 4   #IPENCIL 5   #VAR 6   #OPT_DECOMP 7                 @                            2     PRESENT           
@ @@                              3                     
D @@                             4                      
   @                              5                     
D @@                             6                                  &                   &                   &                                                     
 @@                              7     è             #DECOMP_INFO 	                                                          u #WRITE_SCALAR_REAL 8   #WRITE_SCALAR_COMPLEX =   #WRITE_SCALAR_INTEGER B   #         @     @X                            8                   #FH 9   #DISP :   #N ;   #VAR <             
@ @@                              9                     
D @@                             :                      
   @                              ;                    
@ @@                             <                    
 !   p          5  p        r ;       5  p        r ;                     #         @     @X                            =                   #FH >   #DISP ?   #N @   #VAR A             
@ @@                              >                     
D @@                             ?                      
   @                              @                    
@ @@                             A                     "   p          5  p        r @       5  p        r @                     #         @     @X                            B                   #FH C   #DISP D   #N E   #VAR F             
@ @@                              C                     
D @@                             D                      
   @                              E                    
@ @@                              F                     #   p          5  p        r E       5  p        r E                                                                            u #READ_SCALAR_REAL G   #READ_SCALAR_COMPLEX L   #READ_SCALAR_INTEGER Q   #         @     @X                            G                   #FH H   #DISP I   #N J   #VAR K             
@ @@                              H                     
D @@                             I                      
@ @@                              J                    
D @@                             K                    
 $    p          5  p        r J       5  p        r J                     #         @     @X                            L                   #FH M   #DISP N   #N O   #VAR P             
@ @@                              M                     
D @@                             N                      
@ @@                              O                    
D @@                             P                     %    p          5  p        r O       5  p        r O                     #         @     @X                            Q                   #FH R   #DISP S   #N T   #VAR U             
@ @@                              R                     
D @@                             S                      
@ @@                              T                    
D @@                              U                     &    p          5  p        r T       5  p        r T                                                                            u #WRITE_PLANE_3D_REAL V   #WRITE_PLANE_3D_COMPLEX ^   #         @     @X                            V                  #WRITE_PLANE_3D_REAL%PRESENT W   #IPENCIL X   #VAR Y   #IPLANE Z   #N [   #FILENAME \   #OPT_DECOMP ]                 @                            W     PRESENT           
   @                              X                     
@ @@                             Y                   
 '             &                   &                   &                                                     
   @                              Z                     
   @                              [                     
@ @@                             \                    1           
 @@                              ]     è             #DECOMP_INFO 	   #         @     @X                            ^                  #WRITE_PLANE_3D_COMPLEX%PRESENT _   #IPENCIL `   #VAR a   #IPLANE b   #N c   #FILENAME d   #OPT_DECOMP e                 @                            _     PRESENT           
   @                              `                     
@ @@                             a                    .             &                   &                   &                                                     
   @                              b                     
   @                              c                     
@ @@                             d                    1           
 @@                              e     è             #DECOMP_INFO 	                                                          u #WRITE_EVERY_REAL f   #WRITE_EVERY_COMPLEX o   #         @     @X                            f                  #WRITE_EVERY_REAL%MOD g   #IPENCIL h   #VAR i   #ISKIP j   #JSKIP k   #KSKIP l   #FILENAME m   #FROM1 n                 @                            g     MOD           
   @                              h                     
   @                             i                   
 5             &                   &                   &                                                     
   @                              j                     
   @                              k                     
   @                              l                     
@ @@                             m                    1           
   @                              n           #         @     @X                            o                  #WRITE_EVERY_COMPLEX%MOD p   #IPENCIL q   #VAR r   #ISKIP s   #JSKIP t   #KSKIP u   #FILENAME v   #FROM1 w                 @                            p     MOD           
   @                              q                     
   @                             r                    E             &                   &                   &                                                     
   @                              s                     
   @                              t                     
   @                              u                     
@ @@                             v                    1           
   @                              w                             @               @           	     'è                   #XST x   #XEN y   #XSZ z   #YST {   #YEN |   #YSZ }   #ZST ~   #ZEN    #ZSZ    #X1DIST    #Y1DIST    #Y2DIST    #Z2DIST    #X1CNTS    #Y1CNTS    #Y2CNTS    #Z2CNTS    #X1DISP    #Y1DISP    #Y2DISP    #Z2DISP    #X1COUNT    #Y1COUNT    #Y2COUNT    #Z2COUNT    #EVEN                                                x                                p          p            p                                                                      y                               p          p            p                                                                      z                               p          p            p                                                                      {            $                   p          p            p                                                                      |            0                   p          p            p                                                                      }            <                   p          p            p                                                                      ~            H                   p          p            p                                                                                  T                   p          p            p                                                                                  `              	     p          p            p                                                                                p              
               &                                                                                                 ¸                             &                                                                                                                              &                                                                                                 H                            &                                                                                                                             &                                                                                                 Ø                            &                                                                                                                              &                                                                                                 h                            &                                                                                                 °                            &                                                                                                 ø                            &                                                                                                 @                            &                                                                                                                             &                                                                                            Ð                                                             Ô                                                             Ø                                                             Ü                                                             à                                                                                                                                                                                                                      17#         @                                                    #DECOMP                @                                   è              #DECOMP_INFO 	                                                                                                       22                                                                                                                                                                                     p          p            p                                                                                               p          p            p                                                                                               p          p            p                                                                                               p          p            p                                                                                               p          p            p                                                                                               p          p            p                                       fn#fn "   ¼   ×   b   uapp(DECOMP_2D_IO      @   J  DECOMP_2D    Ó  @   J  MPI (     k       gen@DECOMP_2D_WRITE_ONE    ~        WRITE_ONE_REAL '     @      WRITE_ONE_REAL%PRESENT '   V  @   a   WRITE_ONE_REAL%IPENCIL #     ¼   a   WRITE_ONE_REAL%VAR (   R  L   a   WRITE_ONE_REAL%FILENAME *     Y   a   WRITE_ONE_REAL%OPT_DECOMP "   ÷        WRITE_ONE_COMPLEX *     @      WRITE_ONE_COMPLEX%PRESENT *   Ò  @   a   WRITE_ONE_COMPLEX%IPENCIL &     ¼   a   WRITE_ONE_COMPLEX%VAR +   Î  L   a   WRITE_ONE_COMPLEX%FILENAME -     Y   a   WRITE_ONE_COMPLEX%OPT_DECOMP '   s  i       gen@DECOMP_2D_READ_ONE    Ü        READ_ONE_REAL &   s  @      READ_ONE_REAL%PRESENT &   ³  @   a   READ_ONE_REAL%IPENCIL "   ó  ¼   a   READ_ONE_REAL%VAR '   ¯	  L   a   READ_ONE_REAL%FILENAME )   û	  Y   a   READ_ONE_REAL%OPT_DECOMP !   T
        READ_ONE_COMPLEX )   î
  @      READ_ONE_COMPLEX%PRESENT )   .  @   a   READ_ONE_COMPLEX%IPENCIL %   n  ¼   a   READ_ONE_COMPLEX%VAR *   *  L   a   READ_ONE_COMPLEX%FILENAME ,   v  Y   a   READ_ONE_COMPLEX%OPT_DECOMP (   Ï  k       gen@DECOMP_2D_WRITE_VAR    :        WRITE_VAR_REAL '   Ö  @      WRITE_VAR_REAL%PRESENT "     @   a   WRITE_VAR_REAL%FH $   V  @   a   WRITE_VAR_REAL%DISP '     @   a   WRITE_VAR_REAL%IPENCIL #   Ö  ¼   a   WRITE_VAR_REAL%VAR *     Y   a   WRITE_VAR_REAL%OPT_DECOMP "   ë        WRITE_VAR_COMPLEX *     @      WRITE_VAR_COMPLEX%PRESENT %   Ê  @   a   WRITE_VAR_COMPLEX%FH '   
  @   a   WRITE_VAR_COMPLEX%DISP *   J  @   a   WRITE_VAR_COMPLEX%IPENCIL &     ¼   a   WRITE_VAR_COMPLEX%VAR -   F  Y   a   WRITE_VAR_COMPLEX%OPT_DECOMP '     i       gen@DECOMP_2D_READ_VAR            READ_VAR_REAL &   £  @      READ_VAR_REAL%PRESENT !   ã  @   a   READ_VAR_REAL%FH #   #  @   a   READ_VAR_REAL%DISP &   c  @   a   READ_VAR_REAL%IPENCIL "   £  ¼   a   READ_VAR_REAL%VAR )   _  Y   a   READ_VAR_REAL%OPT_DECOMP !   ¸        READ_VAR_COMPLEX )   V  @      READ_VAR_COMPLEX%PRESENT $     @   a   READ_VAR_COMPLEX%FH &   Ö  @   a   READ_VAR_COMPLEX%DISP )     @   a   READ_VAR_COMPLEX%IPENCIL %   V  ¼   a   READ_VAR_COMPLEX%VAR ,     Y   a   READ_VAR_COMPLEX%OPT_DECOMP +   k         gen@DECOMP_2D_WRITE_SCALAR "   ö  j      WRITE_SCALAR_REAL %   `  @   a   WRITE_SCALAR_REAL%FH '      @   a   WRITE_SCALAR_REAL%DISP $   à  @   a   WRITE_SCALAR_REAL%N &      ´   a   WRITE_SCALAR_REAL%VAR %   Ô  j      WRITE_SCALAR_COMPLEX (   >  @   a   WRITE_SCALAR_COMPLEX%FH *   ~  @   a   WRITE_SCALAR_COMPLEX%DISP '   ¾  @   a   WRITE_SCALAR_COMPLEX%N )   þ  ´   a   WRITE_SCALAR_COMPLEX%VAR %   ²  j      WRITE_SCALAR_INTEGER (     @   a   WRITE_SCALAR_INTEGER%FH *   \  @   a   WRITE_SCALAR_INTEGER%DISP '     @   a   WRITE_SCALAR_INTEGER%N )   Ü  ´   a   WRITE_SCALAR_INTEGER%VAR *            gen@DECOMP_2D_READ_SCALAR !     j      READ_SCALAR_REAL $     @   a   READ_SCALAR_REAL%FH &   Â  @   a   READ_SCALAR_REAL%DISP #      @   a   READ_SCALAR_REAL%N %   B   ´   a   READ_SCALAR_REAL%VAR $   ö   j      READ_SCALAR_COMPLEX '   `!  @   a   READ_SCALAR_COMPLEX%FH )    !  @   a   READ_SCALAR_COMPLEX%DISP &   à!  @   a   READ_SCALAR_COMPLEX%N (    "  ´   a   READ_SCALAR_COMPLEX%VAR $   Ô"  j      READ_SCALAR_INTEGER '   >#  @   a   READ_SCALAR_INTEGER%FH )   ~#  @   a   READ_SCALAR_INTEGER%DISP &   ¾#  @   a   READ_SCALAR_INTEGER%N (   þ#  ´   a   READ_SCALAR_INTEGER%VAR *   ²$  u       gen@DECOMP_2D_WRITE_PLANE $   '%  °      WRITE_PLANE_3D_REAL ,   ×%  @      WRITE_PLANE_3D_REAL%PRESENT ,   &  @   a   WRITE_PLANE_3D_REAL%IPENCIL (   W&  ¼   a   WRITE_PLANE_3D_REAL%VAR +   '  @   a   WRITE_PLANE_3D_REAL%IPLANE &   S'  @   a   WRITE_PLANE_3D_REAL%N -   '  L   a   WRITE_PLANE_3D_REAL%FILENAME /   ß'  Y   a   WRITE_PLANE_3D_REAL%OPT_DECOMP '   8(  ³      WRITE_PLANE_3D_COMPLEX /   ë(  @      WRITE_PLANE_3D_COMPLEX%PRESENT /   +)  @   a   WRITE_PLANE_3D_COMPLEX%IPENCIL +   k)  ¼   a   WRITE_PLANE_3D_COMPLEX%VAR .   '*  @   a   WRITE_PLANE_3D_COMPLEX%IPLANE )   g*  @   a   WRITE_PLANE_3D_COMPLEX%N 0   §*  L   a   WRITE_PLANE_3D_COMPLEX%FILENAME 2   ó*  Y   a   WRITE_PLANE_3D_COMPLEX%OPT_DECOMP *   L+  o       gen@DECOMP_2D_WRITE_EVERY !   »+  ²      WRITE_EVERY_REAL %   m,  <      WRITE_EVERY_REAL%MOD )   ©,  @   a   WRITE_EVERY_REAL%IPENCIL %   é,  ¼   a   WRITE_EVERY_REAL%VAR '   ¥-  @   a   WRITE_EVERY_REAL%ISKIP '   å-  @   a   WRITE_EVERY_REAL%JSKIP '   %.  @   a   WRITE_EVERY_REAL%KSKIP *   e.  L   a   WRITE_EVERY_REAL%FILENAME '   ±.  @   a   WRITE_EVERY_REAL%FROM1 $   ñ.  µ      WRITE_EVERY_COMPLEX (   ¦/  <      WRITE_EVERY_COMPLEX%MOD ,   â/  @   a   WRITE_EVERY_COMPLEX%IPENCIL (   "0  ¼   a   WRITE_EVERY_COMPLEX%VAR *   Þ0  @   a   WRITE_EVERY_COMPLEX%ISKIP *   1  @   a   WRITE_EVERY_COMPLEX%JSKIP *   ^1  @   a   WRITE_EVERY_COMPLEX%KSKIP -   1  L   a   WRITE_EVERY_COMPLEX%FILENAME *   ê1  @   a   WRITE_EVERY_COMPLEX%FROM1 &   *2  o      DECOMP_INFO+DECOMP_2D *   3     a   DECOMP_INFO%XST+DECOMP_2D *   54     a   DECOMP_INFO%XEN+DECOMP_2D *   Ñ4     a   DECOMP_INFO%XSZ+DECOMP_2D *   m5     a   DECOMP_INFO%YST+DECOMP_2D *   	6     a   DECOMP_INFO%YEN+DECOMP_2D *   ¥6     a   DECOMP_INFO%YSZ+DECOMP_2D *   A7     a   DECOMP_INFO%ZST+DECOMP_2D *   Ý7     a   DECOMP_INFO%ZEN+DECOMP_2D *   y8     a   DECOMP_INFO%ZSZ+DECOMP_2D -   9     a   DECOMP_INFO%X1DIST+DECOMP_2D -   ©9     a   DECOMP_INFO%Y1DIST+DECOMP_2D -   =:     a   DECOMP_INFO%Y2DIST+DECOMP_2D -   Ñ:     a   DECOMP_INFO%Z2DIST+DECOMP_2D -   e;     a   DECOMP_INFO%X1CNTS+DECOMP_2D -   ù;     a   DECOMP_INFO%Y1CNTS+DECOMP_2D -   <     a   DECOMP_INFO%Y2CNTS+DECOMP_2D -   !=     a   DECOMP_INFO%Z2CNTS+DECOMP_2D -   µ=     a   DECOMP_INFO%X1DISP+DECOMP_2D -   I>     a   DECOMP_INFO%Y1DISP+DECOMP_2D -   Ý>     a   DECOMP_INFO%Y2DISP+DECOMP_2D -   q?     a   DECOMP_INFO%Z2DISP+DECOMP_2D .   @  H   a   DECOMP_INFO%X1COUNT+DECOMP_2D .   M@  H   a   DECOMP_INFO%Y1COUNT+DECOMP_2D .   @  H   a   DECOMP_INFO%Y2COUNT+DECOMP_2D .   Ý@  H   a   DECOMP_INFO%Z2COUNT+DECOMP_2D +   %A  H   a   DECOMP_INFO%EVEN+DECOMP_2D !   mA  p       MYTYPE+DECOMP_2D $   ÝA  r       REAL_TYPE+DECOMP_2D *   OB  T       GET_DECOMP_INFO+DECOMP_2D 1   £B  Y   e   GET_DECOMP_INFO%DECOMP+DECOMP_2D '   üB  r       COMPLEX_TYPE+DECOMP_2D '   nC  @       MYTYPE_BYTES+DECOMP_2D     ®C  @       NRANK+DECOMP_2D !   îC         XSTART+DECOMP_2D    D         XEND+DECOMP_2D !   E         YSTART+DECOMP_2D    ªE         YEND+DECOMP_2D !   >F         ZSTART+DECOMP_2D    ÒF         ZEND+DECOMP_2D 