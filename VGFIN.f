C***********************************************************************
      program VGFIN
C***********************************************************************
C**** Program to set up a particle dipole array for the VGF code.      *
C****   A particle (sphere, oblate or prolate spheroid, cubical box or *
C****   spheroidal shell) is split into a cubic array of dipole cells  *
C****   (cubic cells with a dipole at the center).  A coarse array of  *
C****   NSID^3 dipole cells is defined using a cubic lattice. The      *
C****   particle is carved out of a NSID^3 cube of dipole cells. The   *
C****   coarse array cells are divided into NLSID fine array cells.    *
C****   Each of the fine array cells is tested to see if it is within  *
C****   the particle geometry. If so, it is counted.  The coarse array *
C****   cell is weighted  in volume (by adjusting the cell length on a *
C****   side (D) according to how many of its fine array cells are in  *
C****   the particle.                                                  *
C****                                                                  *
C**** Units:                                                           *
C****   Angles are input in degrees and converted to radians.          *
C****   Lengths are all relative, that is, if you input a wavelength   *
C****   of 10 um then all other lengths must be in um. You may use     *
C****   m or cm or furlongs if you wish, as long as all lengths are    *
C****   in the same units.                                             *
C****                                                                  *
C**** Time dependence:                                                 *
C****   The VGF code uses the time dependence exp(-iwt). This code     *
C****   (VGFIN) receives the index of refraction from the user and     *
C****   modifies the sign of the imaginary part to fit the VGF time    *
C****   time dependence convention (i.e. the imaginary part should be  *
C****   positive).                                                     *
C****                                                                  *
C**** Geometry:                                                        *
C****   The particle is aligned with its symmetry axis along           *
C****   the Z-axis.  The incident plane wave may be rotated            *
C****   around the particle by the angles alpha, beta, and gamma       *
C****   (corresponding to the Euler angles as given in Arfkin,         *
C****   Mathematical Methods for Physicists, 2nd Ed., Academic Press,  *
C****   NY, 1970). Gamma is effectively the polarization angle around  *
C****   the incident (k-hat) direction and is input as such.           *
C****                                                                  *
C**** Author: Original code by R. T. Lines                             *
C****                                                                  *
C**** extensive modification by Lines 8 APR 96                         *
C****  Last Modification by Lines 17 JAN 98                            *
C***********************************************************************
C**** Variable Definitions                                          ****
      implicit none
      include 'nmax.inc'          
C *** Indices     
      integer G, H, I, IX, IY, IZ
C *** Counters   
      integer  NSID, NLSID, NCUB, INCNT, NUSE, IHP, IC
C *** Reals    
      real W, RAD, alpha, beta, gamma, psi, MR, MI, AR
      real R(NMAX, 3), RF(3), D(NMAX), ER, EI
      real PI, TD, DRF, LD, FRAC, SDN, SD
      real x,y,z     
      real temp
C *** Complex    
      complex M, EPS, C1
C *** Logical (function)    
      logical ISINSIDE
C *** Character    
      character*80 OUTFILE, COMMENTS
C *** Parameters    
      parameter (PI = 3.141592654)
      C1=cmplx(1.0,0.0)                               
      SDN=0.0
C**** Start the program by asking the user for data                 ****
      write(6,*) 'Enter (0) for sphereoidal particles, (1) for cubes'
      read(5,*) IC
      write(6, *) 'Enter wavelength:'
      read(5, *) W
      write(6, *) 'Enter number of dipoles along major axes:'   
      read(5, *) NSID
      write(6, *) 'Enter number of fine structure divisions'
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
      write(6,*) 'Enter choice of solid sphere (0) or hollow shell (1)'
      read(5,*) IHP
      if (IHP.eq.1) then
        write(6,*) 'Enter number of skin depths for thickness of shell'
        read(5,*) SDN    
      end if
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
C *** The number if cells within the particle geometry is counted in
C *** Nuse. This will be equal to Nc from the text.
C
      NUSE = 0           
C
C *** Find the center positions of the fine array cells. Test to  see **
C *   if they are in the particle, if so calculate the volume weight   *
C *   if none of the coarse cell's fine array cells are within the     *
C *   particle geometry, then drop the coarse array cell.              *
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
                  if (ISINSIDE(RAD,RF,AR,SD,SDN,IHP,IC)) INCNT=INCNT+1
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
C*-------------------------------------------------------------------
      logical function ISINSIDE(RAD, RF, AR, SD, SDN, IHP,IC)
C*-------------------------------------------------------------------
C***********************************************************************
C*    Function to determine if a fine array cell is inside the         *
C*      particle. The fine array cell is considered inside if its      *
C*      center point is within the particle.  Only spheroids and       *
C*      cubical boxes are considered here. This version is modified to *
C*      exclude cells within a number of skin depths to allow creation *
C*      of spherical shells if the flag IHP is set.                    *
C***********************************************************************
C**** VAriables                                                     ****
      real RAD, RF(3),AR, SD, SDN
      real RAD2,ARAD2,TRAD
      real SDRAD,SDRAD2,SDARAD2,SDTRAD  
      integer IHP, IC  
C *** If the particle is a cubical box, then just return             ***
      If (IC .eq. 1) then 
         if ((RF(3).le.RAD/AR).and(RF(3).ge.-RAD/AR)) then
            ISINSIDE = .TRUE.
         end if  
         return
      end if
C *** Calculate a semi-major axes that is several (SDN) skin depths  ***
C *     less than the given semi-major axes.                           *
      SDRAD=RAD-SDN*SD
C *** Use the spheroid equation to see if the point falls inside the ***
C *      spheroidal particle. Set up the equation here.                *
      RAD2=RAD*RAD
      ARAD2=RAD2/(AR*AR)
      TRAD = RF(1)*RF(1)/ARAD2                  ! spheroid equation
      TRAD = TRAD + RF(2)*RF(2)/ARAD2
      TRAD = TRAD + RF(3)*RF(3)/RAD2
      TRAD = (TRAD)**0.5
C *** Now if the particle radius is larger than (SDN) skin depths,   ***
C *       use the spheroidal equation to define a smaller spheroid     *
C *       inside the particle and test to see if the point is outside  *
C *       this spheroid.                                               *
      if ((SDRAD.gt.0).and.(IHP.eq.1)) then
         SDRAD2=SDRAD*SDRAD
         SDARAD2=SDRAD2/(AR*AR)
         SDTRAD = RF(1)*RF(1)/SDARAD2           ! spheroid equation
         SDTRAD = SDTRAD + RF(2)*RF(2)/SDARAD2
         SDTRAD = SDTRAD + RF(3)*RF(3)/SDRAD2
         SDTRAD = (SDTRAD)**0.5
      else
         SDTRAD=1.0
      end if
C *** See if it is between the two spheroids, watch out for roundoff.***
      ISINSIDE = ((TRAD.LE.1 + 1.0e-6).AND.(SDTRAD.GE.1 - 1.0e-6))
      return
      end
C***********************************************************************
C******************************************************************************
C*    Include file for maximum array size parameter.
C******************************************************************************
      integer NMAX
      parameter (NMAX = 3000)       
      integer KNMAX
      parameter (KNMAX = 200)  
C******************************************************************************