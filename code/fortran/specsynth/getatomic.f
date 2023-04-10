c------------------------------------------------------------------------------
c
      SUBROUTINE          getatomic
c
c     grab the atomic data file; loop over ions and stuff the atomic
c     data arrays
c
c..............................................................................
c
      include           'specsynth.h'
      include           'const.dek'
      integer           i,linei,natoms
      double precision  w(maxlines),f(maxlines),gamma(maxlines)
      character*80      atom_list,mf_path,def_file,ion_str(maxlines)



      def_file  = 
     & '/fs1/project/cgm_world/code/fortran/input_files/atoms.dat' 



c     the constants assume that the column will be in 1.e13 units and
c     the wavlengths are in angstrom units

      natoms = 0

      OPEN(unit=1,file=def_file,err=1999,status='old')

 01   do 11 i=1,1000
       read(1,*,end=18) ion_str(i),w(i),f(i),gamma(i)
       natoms = natoms + 1
 11   continue
 18   CLOSE(unit=1)


       do 13 linei=1,lines
         do 14 i=1,natoms
          if (ion_name(linei).eq.ion_str(i)) then
           lambda0(linei) = w(i)
           con1(linei)    = 1.0d+5*f(i)*w(i)**2 * sqrt(pi)*e*e/(me*c*c)
           con2(linei)    = 1.0e-8*gamma(i)*w(i)**2 / (4.0d0*pi*c)
           GOTO 13
          end if
 14     continue
 13    continue


      return

c     error opening local parameters, use default and warn


 1999 write(STDOUT,*) ' NOT FOUND- atomic params : atoms.dat'
      STOP ' *** SPECSYNTH(getatomic): terminated ***'
c
      end

c.............................................................................
c

      SUBROUTINE getHIseries
c
c
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'
      integer           i,j,k

      character*80      HIseries(32)


      data HIseries / 'HI1026', 'HI972' , 'HI950' , 'HI938' , 'HI931' , 
     @                'HI926' , 'HI923' , 'HI921' , 'HI919' , 'HI918' , 
     @                'HI917' , 'HI916a', 'HI916b', 'HI915a', 'HI915b',
     @                'HI915c', 'HI914a', 'HI914b', 'HI914c', 'HI914d',    
     @                'HI913a', 'HI913b', 'HI913c', 'HI913d', 'HI913e',    
     @                'HI913f', 'HI913g', 'HI913g', 'HI913h', 'HI913i',    
     @                'LLS'   , 'Lya'   
     @              /

      j = lines
      k = lines

      DO 11 i=1,32
        ion_name(k) = HIseries(i)
        zline(k)    = zline(j)
        nline(k)    = nline(j)
        bline(k)    = bline(j)
        k           = k + 1
 11   CONTINUE

      RETURN
      END
c
c.............................................................................
