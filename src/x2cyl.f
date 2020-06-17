C***********************************************************************
      program CHIIN
C***********************************************************************
C**** Program to set up a particle dipole array for the IBM code.      *
C****   A particle (sphere, oblate or prolate spheroid) is split       *
C****   into a cubic array of dipole cells (cubic cells with a dipole  *
C****   at the center).  A course array of NSID^3 dipole cells is      *
C****   defined using a cubic lattice. The particle is, in effect,     *
C****   carved out of a NSID^3 cube of dipole cells. The course array  *
C****   cells are divided into NLSID fine array cells.  Each of the    *
C****   fine array cells is tested to see if it is within the particle *
C****   geometry. If so, it is counted.  The course array cell is      *
C****   weighted  in volume (buy adjusting the cell length on a side,  *
C****   (D) according to how many of it's fine array cells are in the  *
C****   particle.                                                      *
C
C      HOLE PUNTCHING IS DISABLED
C
C****                                                                  *
C**** Units:                                                           *
C****   Angles are input in degrees and converted to radians.          *
C****   Lengths are all relative, that is, if you input a wavelength   *
C****   of 10 um then all other lengths must be in um. You may use     *
C****   m or cm or furlongs if you wish as long as all lengths are     *
C****   in the same units.                                             *
C****                                                                  *
C**** Time dependence:                                                 *
C****   The IBM code uses the time dependence exp(-iwt). This code     *
C****   (IBMIN) receives the index of refraction from the user and     *
C****   modifies the sign of the imaginary part to fit the IBM time    *
C****   time dependence convention (i.e. the imaginary part should be  *
C****   positive).                                                     *
C****                                                                  *
C**** Geometry:                                                        *
C****   The particle is aligned with it's symmetry axis along          *
C****   the Z-axis.  The incident plane wave may be rotated            *
C****   around the particle by the angles alpha, beta, and gamma       *
C****   (corresponding to the Euler angles as given in Arfkin,         *
C****   Mathematical Methods for Physicists, 2nd Ed., Academic Press,  *
C****   NY, 1970). Gamma is effectively the polarization angle around  *
C****   the incident (k-hat) direction and is input as such.           *
C****                                                                  *
C**** Author: Original code by Lines                                   *
C****                                                                  *
C**** extensive modification by Lines 8 APR 96                         *
C****  Last Modification by Lines 26 DEC 96                            * 
C***********************************************************************
C**** Variable Definitions                                          ****
      implicit none
      include 'nmax.inc'           
C *** Indices      
      integer G, H, I, IX, IY, IZ
C *** Counters    
      integer  NSID, NLSID, NCUB, INCNT, NUSE
C *** Reals     
      real W, RAD, alpha, beta, gamma, psi, MR, MI, AR
      real R(NMAX, 3), RF(3), D(NMAX), ER, EI
      real PI, TD, DRF, LD, FRAC, SD
      real x,y,z      
      real temp
C *** Complex     
      complex M, EPS, C1
C *** Logical (function)     
      logical ISINSIDE2
C *** Character     
      character*80 OUTFILE, COMMENTS
C *** Parameters     
      parameter (PI = 3.141592654)
      C1=cmplx(1.0,0.0)
C**** Start the program by asking the user for data                 ****      
      write(6, *) 'Enter wavelength:'
      read(5, *) W
      write(6, *) 'Enter number of dipoles along major axes:'    
      read(5, *) NSID
      write(6, *) 'Enter number of fine structure devisions' 
      write(6, *) 'per dipole cell side:'
      read(5, *) NLSID 
      write(6, *) 'Enter particle symmetry semi-axes (z-axes radius):'
      read(5, *) RAD
      write(6, *) 'Enter incident theta, phi, gamma (Degrees):'      
      read(5, *) alpha, beta, gamma
      write(6, *) 'Enter polarization angle (Degrees):'      
      read(5, *) psi
      write(6, *) 'Enter real & imag parts of index of refraction:'
      read(5, *) MR, MI
      write(6, *) 'Enter aspect ratio (major axis/minor axis)'
      read(5,*) AR
      write(6, *) 'Enter file comments:'
      read(5, '(1a80)') COMMENTS
      write(6, *) 'Enter output filename:'      
      read(5, '(1a80)') OUTFILE
      write(6, *)        
C**** Start calculations                                            ****
C *** Calculate skin depth
      if (MI.ne.0.0) then
         SD=W/(2*pi*MI)
      else 
         SD=0.0
      end if
C *** Get dipole coarse array cell length on a side i                ***
      TD = 2.0*max(RAD,RAD/AR)/NSID  
C *** Make sure the imaginary index of refraction is the right sign  ***
C *     We will use exp(ikr-iwt) in ADMF so MI positive                *
      if (MI.LT.0.0) MI = -MI
      M = cmplx(MR, MI)        
C *** Half a cell dimension ***
      DRF = 0.5*(float(NSID) - 1.0)
C *** Build coarse array... (first let user know what is happening)  ***     
      write(6, *) 'Building coarse array...'
      NCUB = 0               ! Number of total cells (N-cubed)
C *** Fine array cell size (little D)                                ***
      LD = TD/NLSID
C *** 1 over Volume of the fine array cell (used later)              ***      
      FRAC = (1.0/float(NLSID))**3
C
      write(6, *) 'Calculating fine structure...'
C
      NUSE = 0            ! Number of cells within particle geometry
C
C *** Find the center positions of the fine array cells. Test to  see **
C *   if they are in the particle, if so calculate the volume weight   *
C *   if not then drop the coarse array cell.                          *
      do G = 0, NSID-1
        do H = 0, NSID-1
          do I = 0, NSID-1            
            x = (float(G) - DRF)*TD
            y = (float(H) - DRF)*TD
            z = (float(I) - DRF)*TD
            INCNT = 0    
            do IX = 0, NLSID-1
              do IY = 0, NLSID-1
                do IZ = 0, NLSID-1
                  RF(1) = x + LD*(0.5*(1.0 - NLSID) + IX)      
                  RF(2) = y + LD*(0.5*(1.0 - NLSID) + IY)
                  RF(3) = z + LD*(0.5*(1.0 - NLSID) + IZ)
                  if (ISINSIDE2(RAD, RF, AR, SD)) INCNT=INCNT+1
                end do   !IZ
              end do     !IY
            end do       !IX
            if (INCNT.GT.0) then                ! drop unfilled cells      
              NUSE = NUSE + 1       
              if (MOD(NUSE, 1000).EQ.0) write(6, *) NUSE
              R(NUSE, 1) = x 
              R(NUSE, 2) = y
              R(NUSE, 3) = z
              D(NUSE) = TD * (INCNT*FRAC)**(1.0/3.0)   !weight cells
            endif
          end do
        end do
      end do                    
C**** DONE! Write out the cubic array and data to file              **** 
      write(6, *) 'Writing file...'
30    open(10, file = OUTFILE)
      write(10, *) COMMENTS
      write(10, *) NUSE
      write(10, *) W, alpha, beta, gamma, psi, RAD      
      EPS = M*M         
      ER = real(EPS)
      EI = aimag(EPS)
      write(10,*) ER, EI, TD
      do I = 1, NUSE
        write(10, 100) R(I, 1), R(I, 2), R(I, 3), D(I)
        temp=temp+d(I)**3
      end do
      write(10,*) temp
      write(10,*) 4.*PI*RAD**3/3.
      close(10)
100   format(4(1x, 1e12.4))
      end

C***********************************************************************
C*---------------------------------------------------------------------*
      logical function ISINSIDE2(RAD, RF, AR, SD)
C*---------------------------------------------------------------------*
C***********************************************************************   
C*    Function to determine if a fine array cell is inside the         *
C*      The fine array cell is considered inside if it's center point  *
C*      is within the particle.  Only spheroids are considered here.   *
C*      This version is modified to exclude cells within a number of   *
C*      skin depths.                                                   *
C***********************************************************************   
C**** VAriables                                                     ****
      real RAD, RF(3),AR, SD 
      real RAD2,ARAD2,TRAD
      real SDRAD,SDRAD2,SDARAD2,SDTRAD   
C *** Calculate a semi-major axes that is several (3) skin depths    ***
C *     less than the given semi-major axes.                           *
      SDRAD=RAD-3*SD
C *** Use the spheroid equation to see if the point falls inside the ***
C *      spheroidal particle. Set up the equation here.                *
      RAD2=RAD*RAD
      ARAD2=RAD2/(AR*AR)
      TRAD = RF(1)*RF(1)/ARAD2                  ! spheroid equation 
      TRAD = TRAD + RF(2)*RF(2)/ARAD2
      TRAD = TRAD + RF(3)*RF(3)/RAD2
      TRAD = (TRAD)**0.5
C *** Now if the particle radius is larger than (3) skin depths,     ***
C *       use the spheroidal equation to define a smaller spheroid     *
C *       inside the particle and test to see if the point is outside  *
C *       this spheroid.                                               *
c      if (SDRAD.gt.0) then
c         SDRAD2=SDRAD*SDRAD
c         SDARAD2=SDRAD2/(AR*AR)
c         SDTRAD = RF(1)*RF(1)/SDARAD2           ! spheroid equation 
c         SDTRAD = SDTRAD + RF(2)*RF(2)/SDARAD2
c         SDTRAD = SDTRAD + RF(3)*RF(3)/SDRAD2
c         SDTRAD = (SDTRAD)**0.5
c      else
         SDTRAD=1.0
c      end if
C *** See if it is between the two spheroids, watch out for roundoff.***
      ISINSIDE2 = ((TRAD.LE.1 + 1.0e-6).AND.(SDTRAD.GE.1 - 1.0e-6)) 
      return
      end
C***********************************************************************
