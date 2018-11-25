module global_toolkit_module
!
!  Purpose:
!   this module contains a collection of useful procedures used in several other
!   modules.
!
!   **** PROGRAMMING PRINCIPLES: ****
!   all procedures are PURE.
!   FIX the sizes of dummy args whenever possible.
!     - suffix '2d'/'3d' is added at the end of procedure name to indicate
!       its supported dimension
!     - the sizes can be set by private parameters in the procedure.
!   AVOID support of variable sizes of dummy args
!     - Keep It Simple and Stupid (KISS)
!     - forbid unexpected input arg. sizes, avoid explicit checking
!   AVOID optional arg. KISS
!   if error status and message are needed, make them COMPULSARY INPUTS
!     - do NOT give the option to not check for error status
!   ALWAYS CHECK dummy arg sizes if variable sizing has to be allowed,
!       and flag error when their sizes are not expected
!   ALWAYS CHECK the validity of intent in/inout dummy args if they are expected
!       to be within a certain range of values
!       EXCEPTION: if the input arg is a size parameter of other input args,
!       then its checking can be left to the compilor
!   ALWAYS USE LOCAL COPY OF INTENT(INOUT) and OPTIONAL dummy args. for
!       calculation, only update to the original arg. before successful return
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    15/04/15  B. Y. Chen            Original code
!

implicit none

contains



  pure function glb_dee3d (dee, theta)
  ! Purpose :
  ! to transform d matrix from local to global coordinates
  ! (for composite lamina, 3D, with fibre angle = theta)

    use parameter_module, only : DP, ZERO, PI, HALFCIRC, TWO, FOUR

    ! define private parameters
    ! NST: no. of stress/strain terms, used to set the size of dummy args
    integer, parameter   :: NST = 6

    real(DP), intent(in) :: dee(NST,NST), theta
    real(DP)             :: glb_dee3d(NST,NST)

    ! local variables
    real(DP) :: c, s

    glb_dee3d = ZERO
    c = ZERO
    s = ZERO

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)

    ! d matrix in global coords stored first in local array glb_dee3d
    glb_dee3d(1,1) = c*c*c*c*dee(1,1) + TWO * c*c*s*s*(dee(1,2)       &
                 & + TWO * dee(4,4)) + s*s*s*s*dee(2,2)
    glb_dee3d(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - FOUR * dee(4,4))  &
                 & + (s*s*s*s+c*c*c*c)*dee(1,2)
    glb_dee3d(2,1) = glb_dee3d(1,2)
    glb_dee3d(2,2) = s*s*s*s*dee(1,1) + TWO * c*c*s*s*(dee(1,2)       &
                 & + TWO * dee(4,4)) + c*c*c*c*dee(2,2)
    glb_dee3d(1,3) = c*c*dee(1,3) + s*s*dee(2,3)
    glb_dee3d(3,1) = glb_dee3d(1,3)
    glb_dee3d(1,4) = s*c*(c*c*(dee(1,1) - dee(1,2) - TWO * dee(4,4))  &
                 & + s*s*(dee(1,2) - dee(2,2) + TWO * dee(4,4)))
    glb_dee3d(4,1) = glb_dee3d(1,4)
    glb_dee3d(2,3) = s*s*dee(1,3) + c*c*dee(2,3)
    glb_dee3d(3,2) = glb_dee3d(2,3)
    glb_dee3d(2,4) = s*c*(s*s*(dee(1,1) - dee(1,2) - TWO * dee(4,4))  &
                 & + c*c*(dee(1,2) - dee(2,2) + TWO * dee(4,4)))
    glb_dee3d(4,2) = glb_dee3d(2,4)
    glb_dee3d(3,3) = dee(3,3)
    glb_dee3d(3,4) = c*s*(dee(1,3) - dee(2,3))
    glb_dee3d(4,3) = glb_dee3d(3,4)
    glb_dee3d(5,5) = c*c*dee(5,5) + s*s*dee(6,6)
    glb_dee3d(5,6) = c*s*(dee(6,6) - dee(5,5))
    glb_dee3d(6,5) = glb_dee3d(5,6)
    glb_dee3d(6,6) = s*s*dee(5,5) + c*c*dee(6,6)
    glb_dee3d(4,4) = s*s*c*c*(dee(1,1) - TWO * dee(1,2) + dee(2,2))   &
                 & + (s*s - c*c)*(s*s - c*c)*dee(4,4)

  end function glb_dee3d




  pure function determinant2d (jacob)
  ! Purpose:
  ! returns the determinant of a 2D jacobian matrix

    use parameter_module, only : DP, ZERO

    ! define private parameters
    ! NDIM: no. of dimensions, used to set the size of dummy arg
    integer, parameter   :: NDIM = 2

    real(DP), intent(in) ::   jacob(NDIM,NDIM)
    real(DP)             ::   determinant2d

    determinant2d = ZERO

    determinant2d = jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1)

  end function determinant2d




  pure function determinant3d (jacob)
  ! Purpose:
  ! returns the determinant of a 3D jacobian matrix

    use parameter_module, only : DP, ZERO

    ! define private parameters
    ! NDIM: no. of dimensions, used to set the size of dummy arg
    integer, parameter   :: NDIM = 3

    real(DP), intent(in) ::   jacob(NDIM,NDIM)
    real(DP)             ::   determinant3d

    determinant3d = ZERO

    determinant3d = jacob(1,1) * &
                & (jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3))
    determinant3d = determinant3d - jacob(1,2) * &
                & (jacob(2,1)*jacob(3,3)-jacob(3,1)*jacob(2,3))
    determinant3d = determinant3d + jacob(1,3) * &
                & (jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2))

  end function determinant3d




  pure subroutine invert_self3d (jacob, istat, emsg, detj)
  ! Purpose:
  ! calculate the inverse of a 3D jacobian matrix and update onto itself
  ! the invertability of jacob must be checked before inverting it

    use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                          & ZERO, SMALLNUM

    ! define private parameters
    ! NDIM: no. of dimensions, used to set the size of dummy arg
    integer, parameter   :: NDIM = 3

    real(DP),                 intent(inout) :: jacob(NDIM,NDIM)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: detj

    ! local variables:
    ! local copy of intent inout arg
    real(DP) :: jacob_lcl(NDIM,NDIM)
    ! inverse of jacob matrix
    real(DP) :: inv_jacob(NDIM,NDIM)
    ! local copy of optional arg.
    real(DP) :: det

    ! initialize intent out and local variables
    istat     = STAT_SUCCESS
    emsg      = ''
    jacob_lcl = ZERO
    inv_jacob = ZERO
    det       = ZERO

    ! assign values to local variables
    jacob_lcl = jacob

    ! get the determinant of jacob matrix; pass local copy of jacob for safety
    if(present(detj)) then
      det = detj
    else
      det = determinant3d(jacob_lcl)
    end if

    ! check to see if the matrix is singular; if so, flag error and exit program
    if(det < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'the jacobian matrix is singular, invert_self3d, &
            & global_toolkit_module'
      return
    end if

    ! calculate the inverse of jacob matrix
    inv_jacob(1,1)=jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3)
    inv_jacob(2,1)=jacob(3,1)*jacob(2,3)-jacob(2,1)*jacob(3,3)
    inv_jacob(3,1)=jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2)
    inv_jacob(1,2)=jacob(3,2)*jacob(1,3)-jacob(1,2)*jacob(3,3)
    inv_jacob(2,2)=jacob(1,1)*jacob(3,3)-jacob(3,1)*jacob(1,3)
    inv_jacob(3,2)=jacob(3,1)*jacob(1,2)-jacob(1,1)*jacob(3,2)
    inv_jacob(1,3)=jacob(1,2)*jacob(2,3)-jacob(2,2)*jacob(1,3)
    inv_jacob(2,3)=jacob(2,1)*jacob(1,3)-jacob(1,1)*jacob(2,3)
    inv_jacob(3,3)=jacob(1,1)*jacob(2,2)-jacob(2,1)*jacob(1,2)

    inv_jacob = inv_jacob/det

    ! update intent(inout) args before successful return
    jacob = inv_jacob

  end subroutine invert_self3d




  pure function beemat3d (gn, nnode)
  ! Purpose:
  ! inputs the gradient of shape functions and
  ! returns the strain-displacement matrix for infinitesimal deformation
  ! order convention of strain terms in the strain vector:
  ! eps1, eps2, eps3, gamma12, gamma13, gamma23

    use parameter_module, only : DP, ZERO

    ! define private parameters, used to set the size of dummy arg
    ! NDIM: no. of dimensions
    ! NST : no. of stress/strain terms
    integer, parameter   :: NDIM = 3, NST = 6

    integer,  intent(in) :: nnode
    real(DP), intent(in) :: gn(nnode, NDIM)
    real(DP)             :: beemat3d(NST, nnode*NDIM)

    ! local variables
    real(DP) :: x, y, z
    integer  :: k, l, m, n

    ! initialize intent out and local variables
    beemat3d = ZERO
    x = ZERO
    y = ZERO
    z = ZERO
    k = 0; l = 0; m = 0; n = 0

    ! nnode is size parameter, let compilor to check it
    ! gn is fixed sized and can take any value, no check

    do m=1, nnode
      n = 3*m
      k = n-1
      l = k-1
      x = gn(m,1)
      y = gn(m,2)
      z = gn(m,3)
      beemat3d(1,l) = x
      beemat3d(4,k) = x
      beemat3d(5,n) = x
      beemat3d(2,k) = y
      beemat3d(4,l) = y
      beemat3d(6,n) = y
      beemat3d(3,n) = z
      beemat3d(5,l) = z
      beemat3d(6,k) = z
    end do

  end function beemat3d




  pure function lcl_strain3d (strain, theta)
  ! Purpose:
  ! transfer strains from global to local coordinate systems
  ! (for 3D composite lamina, ply angle = theta)

    use parameter_module, only : DP, ZERO, PI, HALFCIRC, ONE, TWO

    ! define private parameters, used to set the size of dummy arg
    ! NST : no. of stress/strain terms
    integer, parameter   :: NST = 6

    real(DP), intent(in) :: strain(NST), theta
    real(DP)             :: lcl_strain3d(NST)

    real(DP) :: c, s, Q_matrix(NST,NST)

    lcl_strain3d = ZERO
    c = ZERO
    s = ZERO
    Q_matrix = ZERO

    ! strain and theta can take any value, nothing to check

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)

    ! calculate rotational matrix Q
    Q_matrix(1,1)=c*c
    Q_matrix(1,2)=s*s
    Q_matrix(1,4)=c*s
    Q_matrix(2,1)=s*s
    Q_matrix(2,2)=c*c
    Q_matrix(2,4)=-c*s
    Q_matrix(3,3)=ONE
    Q_matrix(4,1)=-TWO*c*s
    Q_matrix(4,2)=TWO*c*s
    Q_matrix(4,4)=c*c-s*s
    Q_matrix(5,5)=c
    Q_matrix(5,6)=-s
    Q_matrix(6,5)=s
    Q_matrix(6,6)=c

    ! rotate global strain to obtain local strain
    lcl_strain3d = matmul(Q_matrix,strain)

  end function lcl_strain3d




  pure subroutine normalize_vect (a, is_zero_vect, amag)
  ! Purpose:
  ! normalize vector and return its magnitude

  use parameter_module, only : DP, ZERO, SMALLNUM

    real(DP),           intent(inout) :: a(:)
    logical,            intent(out)   :: is_zero_vect
    real(DP), optional, intent(out)   :: amag

    ! local copy of amag
    real(DP) :: mag

    ! initialize intent out and local variables
    is_zero_vect = .false.
    if(present(amag)) amag = ZERO
    mag = ZERO

    ! a can take any value, nothing to check

    ! calculate the length of a
    mag = sqrt(dot_product(a,a))

    if (mag > ZERO) then
    ! a is NOT a zero-length vector;
    ! normalize a, update amag (if present) and return
      is_zero_vect = .false.
      a = a / mag
      if(present(amag)) amag = mag
      return
    else
    ! a is a zero-length vector;
    ! do NOT update a, just return
      is_zero_vect = .true.
      return
    end if

  end subroutine normalize_vect




  pure function cross_product3d (a, b)
  ! Purpose:
  ! function to compute the cross-product of TWO 3D vectors, a * b

  use parameter_module, only : DP

    real(DP), intent(in)  :: a(3), b(3)
    real(DP)              :: cross_product3d(3)

    cross_product3d(1) = a(2) * b(3) - a(3) * b(2)
    cross_product3d(2) = a(3) * b(1) - a(1) * b(3)
    cross_product3d(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product3d




  pure function distance (a, b, n)
  ! Purpose:
  ! function to compute the distance between TWO vectors
  ! their common size n is set as a required input to avoid having to
  ! explicitly check for their size consistency

  use parameter_module, only : DP, ZERO

    integer,  intent(in)  :: n
    real(DP), intent(in)  :: a(n), b(n)
    real(DP)              :: distance
    ! local variable
    real(DP) :: c(n)
    ! initialize local variable
    c = ZERO

    ! intent in variables a and b are fixed-sized and may take any value
    ! intent in variable n must > 0; it is a dummy arg size parameter,
    ! so the checking of n is left to the compilor

    c = a - b

    distance=sqrt(dot_product(c,c))

  end function distance




  pure subroutine edge_crack_cross2d (edge, crack, cross_stat, cross_point)
  ! Purpose:
  ! to checks if TWO line sections cross each other in 2D,
  ! return the cross status and locate the crossing point if they cross

  use parameter_module, only : DP, ZERO, SMALLNUM, &
                        & CROSS_ON_EDGE_ON_CRACK,  CROSS_ON_EDGE_OFF_CRACK,  &
                        & CROSS_OFF_EDGE_ON_CRACK, CROSS_OFF_EDGE_OFF_CRACK, &
                        & EDGE_CRACK_PARALLEL, ENDNODE_CLEARANCE_RATIO

    ! dummy args list:
    ! edge        : nodal coords of two end nodes of edge  line segment
    ! crack       : nodal coords of two end nodes of crack line segment
    ! cross_stat  : status of crossing
    ! cross_point : coords of the cross point
    real(DP), intent(in)  :: edge(2,2)
    real(DP), intent(in)  :: crack(2,2)
    integer,  intent(out) :: cross_stat
    real(DP), intent(out) :: cross_point(2)

    ! local variables
    ! a/b/cLj,      j=1,2 : line equation constants a, b, c of Line j
    ! det                 : determinant of the system equations of two lines
    real(DP) :: aL1, bL1, cL1
    real(DP) :: aL2, bL2, cL2
    real(DP) :: det

    ! initialize intent out and local variables
    ! intent out
    cross_stat  = 0
    cross_point = ZERO
    ! line equation constants and determinant
    aL1 = ZERO; bL1 = ZERO; cL1 = ZERO
    aL2 = ZERO; bL2 = ZERO; cL2 = ZERO
    det = ZERO

    ! intent in dummy args' sizes are fixed, and they may take on any value,
    ! no need to check their validity

    !
    !**** algorithm for finding the intersection of two lines ****
    ! equation of a line: aL * x + bL * y = cL
    ! to find the cross point of line 1 (edge) and line 2 (crack)
    ! is to solve the system equations of two lines:
    ! aL1 * x + bL1 * y = cL1 (equation of line 1)
    ! aL2 * x + bL2 * y = cL2 (equation of line 2)
    ! calculate det = aL1 * bL2 - aL2 * bL1
    ! if det == ZERO, lines are parallel, no intersection possible;
    ! if det /= ZERO, lines are NOT parallel, the cross point coordinates are:
    ! x_crosspoint = (bL2 * cL1 - bL1 * cL2) / det
    ! y_crosspoint = (aL1 * cL2 - aL2 * cL1) / det
    !

    ! line 1 equation (edge)
    call line_equation(edge, aL1, bL1, cL1)

    ! line 2 equation (crack)
    call line_equation(crack, aL2, bL2, cL2)

    ! calculate determinant
    det   = aL1 * bL2 - aL2 * bL1


    !**** MAIN CALCULATIONS ****

    ! select what to do based on the value of det:
    ! if det is NONZERO: intersection is possible
    !   - calculate cross point coordinates
    !   - if it is on edge:
    !       update the cross_stat:
    !         * if it is on crack : cross_stat = CROSS_ON_EDGE_ON_CRACK
    !         * if it is off crack: cross_stat = CROSS_ON_EDGE_OFF_CRACK
    !       adjust cross point coordinates away from edge end nodes,
    !       then exit the program
    !   - if it is off edge:
    !       update the cross_stat:
    !         * if it is on crack : cross_stat = CROSS_OFF_EDGE_ON_CRACK
    !         * if it is off crack: cross_stat = CROSS_OFF_EDGE_OFF_CRACK
    !       ZERO cross point coordinates,
    !       then exit the program
    !
    ! if det is ZERO: edge and crack are parallel, intersection NOT possible
    !   update the cross_stat:
    !     * cross_stat = EDGE_CRACK_PARALLEL
    !   ZERO cross point coordinates,
    !   then exit the program

    det_if: if (abs(det) > ZERO) then
    ! intersection possible

        ! calculate cross point coordinates
        cross_point(1) = (bL2 * cL1 - bL1 * cL2) / det
        cross_point(2) = (aL1 * cL2 - aL2 * cL1) / det

        ! update cross status
        cross_edge_if: if ( in_range(cross_point, edge) ) then
        ! the cross point is on this edge
            if ( in_range(cross_point, crack) ) then
            ! cross point on edge and also on crack
              cross_stat = CROSS_ON_EDGE_ON_CRACK   ! = 2
            else
            ! cross point on edge but not on crack
              cross_stat = CROSS_ON_EDGE_OFF_CRACK  ! = 1
            end if
            ! adjust cross point positions away from edge end nodes
            call adjust_cross_point(cross_point, edge)
        else cross_edge_if
        ! the cross point is out of the range of the edge
            if ( in_range(cross_point, crack) ) then
            ! cross point off edge but on crack
              cross_stat = CROSS_OFF_EDGE_ON_CRACK  ! = -1
            else
            ! cross point off edge and off crack
              cross_stat = CROSS_OFF_EDGE_OFF_CRACK ! = -2
            end if
            ! zero cross point coordinates
            cross_point = ZERO
        end if cross_edge_if

        ! exit program
        return

    else det_if
    ! lines are parallel; no intersection possible

        ! just update output variables and return
        cross_stat  = EDGE_CRACK_PARALLEL           ! = 0
        cross_point = ZERO
        ! exit program
        return

    end if det_if

    !**** END MAIN CALCULATIONS ****


    contains

    ! internal procedures


    pure subroutine line_equation(line_seg, a, b, c)
    ! Purpose:
    ! obtain line equation constants a, b, c with passed-in line segment

      ! dummy args
      real(DP), intent(in)  :: line_seg(2,2)
      real(DP), intent(out) :: a, b, c

      ! local variables:
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L

      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      a   = yN1L - yN2L
      b   = xN2L - xN1L
      c   = yN1L * xN2L - xN1L * yN2L

    end subroutine line_equation


    pure function in_range(point, line_seg)
    ! Purpose:
    ! return true of the point lies within the range of the line segment

      ! dummy args:
      real(DP), intent(in) :: point(2), line_seg(2,2)
      logical              :: in_range

      ! local variables:
      real(DP) :: xp, yp
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L
      real(DP) :: xminL, yminL, xmaxL, ymaxL

      ! coords of point
      xp   = point(1)
      yp   = point(2)
      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      ! range of line_seg
      xminL = min(xN1L, xN2L)
      yminL = min(yN1L, yN2L)
      xmaxL = max(xN1L, xN2L)
      ymaxL = max(yN1L, yN2L)

      ! point is in range of line seg if
      ! xminL <= xp <= xmaxL and
      ! yminL <= yp <= ymaxL
      if (xminL <= xp .and. xp <= xmaxL .and. &
      &   yminL <= yp .and. yp <= ymaxL) then
        in_range = .true.
      else
        in_range = .false.
      end if

    end function in_range


    pure subroutine adjust_cross_point(point, line_seg)
    ! Purpose:
    ! adjust the coords of cross_point to be away from the end nodes of edge

      real(DP), intent(inout) :: point(2)
      real(DP), intent(in)    :: line_seg(2,2)

      ! local variables:
      real(DP) :: xp, yp
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L
      real(DP) :: xminL, yminL, xmaxL, ymaxL

      ! coords of point
      xp   = point(1)
      yp   = point(2)
      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      ! range of line_seg
      xminL = min(xN1L, xN2L)
      yminL = min(yN1L, yN2L)
      xmaxL = max(xN1L, xN2L)
      ymaxL = max(yN1L, yN2L)

      ! adjust point positions away from line seg end nodes
      xp = max(xp, xminL + ENDNODE_CLEARANCE_RATIO * (xmaxL-xminL))
      xp = min(xp, xmaxL - ENDNODE_CLEARANCE_RATIO * (xmaxL-xminL))
      yp = max(yp, yminL + ENDNODE_CLEARANCE_RATIO * (ymaxL-yminL))
      yp = min(yp, ymaxL - ENDNODE_CLEARANCE_RATIO * (ymaxL-yminL))

      ! update point
      point = [xp, yp]

    end subroutine adjust_cross_point


  end subroutine edge_crack_cross2d




  pure subroutine crack_elem_centroid2d (nedge, crack_angle, coords, &
  & nodes_on_edges, istat, emsg, edge_crack_points, crack_edge_indices)
  ! Purpose:
  ! to find TWO cross points between the crack line (passing the centroid) and
  ! the edges of an element, and the two edges crossed.

  use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,  &
                        & ZERO, HALF, HALFCIRC, PI, SMALLNUM, ONE, QUARTER,&
                        & CROSS_ON_EDGE_ON_CRACK,  CROSS_ON_EDGE_OFF_CRACK, &
                        & NTESTCRACKPOINT

    ! list of dummy args:
    ! - nedge              : no. of edges in this element; = no. of nodes
    ! - crack_angle
    ! - coords             : nodal coordinates of this element
    ! - nodes_on_edges(:,i) : indices of the two end nodes of the edge i
    ! - edge_crack_points  : coords of the two crack points on two edges
    ! - crack_edge_indices : indices of the two edges passed by the crack
    integer,  intent(in)  :: nedge
    real(DP), intent(in)  :: crack_angle
    real(DP), intent(in)  :: coords(2, nedge)
    integer,  intent(in)  :: nodes_on_edges(2, nedge)
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    real(DP),                 intent(out) :: edge_crack_points(2,2)
    integer,        optional, intent(out) :: crack_edge_indices(2)

    ! local variables
    ! - centroid          : coordinates of element centroid
    ! - theta             : ply angle in radiant
    ! - crack_unit_vector : unit vector along the crack line
    ! - coords_crack      : coordinates of end nodes of crack line segment
    ! - coords_edge       : coordinates of end nodes of an element edge
    ! - cross_point       : coordinates of the cross point btw crack and edge
    ! - n_crack_edge      : no. of cracked edges
    ! - cross_stat        : status of crossing btw a crack and an edge
    ! - crack_edge_IDs    : local copy of optional arg crack_edge_indices
    real(DP) :: centroid(2)
    real(DP) :: theta, crack_unit_vect(2)
    real(DP) :: coords_crack(2,2), coords_edge(2,2), cross_point(2)
    integer  :: n_crack_edge, cross_stat, crack_edge_IDs(2)
    real(DP) :: projection, maxprojection, minprojection, test_line(2)
    logical  :: is_zero_vect

    integer :: i, j

    ! initialize intent out variables
    istat             = STAT_SUCCESS
    emsg              = ''
    edge_crack_points = ZERO
    if (present(crack_edge_indices)) crack_edge_indices = 0
    ! initialize local variables
    centroid     = ZERO
    theta        = ZERO
    crack_unit_vect = ZERO
    coords_crack = ZERO
    coords_edge  = ZERO
    cross_point  = ZERO
    projection   = ZERO
    maxprojection= ZERO
    minprojection= ZERO
    test_line    = ZERO
    cross_stat     = 0
    n_crack_edge   = 0
    crack_edge_IDs = 0
    is_zero_vect   = .false.
    i=0; j=0

    ! check validity of inputs: nedge, crack_angle, coords, nodes_on_edges
    ! nedge must be > 0, this check is omitted as if it is not satisfied,
    ! compilor should flag error as it sets the size of coords and nodes_on_edges
    ! crack_angle and coords can take any values
    ! nodes_on_edges must be >= 1, this need to be checked
    if ( any(nodes_on_edges < 1) ) then
      istat = STAT_FAILURE
      emsg  = 'end node indices must be >=1, crack_elem_centroid2d, &
      &global_toolkit_module'
      return
    end if

    ! find centroid
    centroid(1) = sum(coords(1,:))/real(nedge,DP)
    centroid(2) = sum(coords(2,:))/real(nedge,DP)

    ! find crack line unit vector
    theta           = crack_angle / HALFCIRC * PI
    crack_unit_vect = [cos(theta), sin(theta)]

    ! set centroid as node 1 of crack line segment
    coords_crack(:,1) = centroid(:)
    ! set node 2 coordinates of crack line segment
    ! to be a unit distance away from node 1 along the crack line
    coords_crack(:,2) = coords_crack(:,1) + crack_unit_vect(:)



    !**** MAIN CALCULATIONS ****

    ! loop over all edges to find TWO cross points with the crack line
    ! theoretically, this should find EXACTLY TWO cross points, as the
    ! crack line passes the element centroid and it must cross TWO edges
    do i = 1, nedge
      ! the two end nodes' coords of edge i
      coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
      coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
      ! zero cross_stat and cross_point for reuse
      cross_stat  = 0
      cross_point = ZERO
      ! check if the edge and crack cross each other
      call edge_crack_cross2d (coords_edge, coords_crack, cross_stat, cross_point)
      ! check cross_stat, update only if cross point is on edge
      if (cross_stat == CROSS_ON_EDGE_ON_CRACK .or. &
      &   cross_stat == CROSS_ON_EDGE_OFF_CRACK) then
          ! increase the no. of cracked edges
          n_crack_edge = n_crack_edge + 1
          ! store the index of this cracked edge
          crack_edge_IDs(n_crack_edge) = i
          ! store the cross point coords of this cracked edge
          edge_crack_points(:, n_crack_edge) = cross_point(:)
       endif
       ! exit the loop when TWO broken edges are already found
       if (n_crack_edge == 2) exit
    end do

    !**** END MAIN CALCULATIONS ****

    ! check if indeed two cracked edges have been found; if not, then the elem
    ! is likely to be poorly shaped, flag an error, clean up and return
    !~if (n_crack_edge /= 2) then
    !~  istat = STAT_FAILURE
    !~  emsg  = 'two cracked edges cannot be found, element is likely to be &
    !~  &poorly shaped, crack_elem_centroid2d, global_toolkit_module'
    !~  ! clean up intent out variable before error exit
    !~  edge_crack_points = ZERO
    !~  return
    !~end if

    select case(n_crack_edge)
    ! if two edges are found, then continue to return
    case (2)
      continue
    ! if only one edge is found, then use an edge midpoint to estimate the other
    case (1)
      ! initialize projection
      projection = ZERO
      ! 1st crack point already found
      coords_crack(:,1) = edge_crack_points(:,1)
      do i = 1, nedge
        ! if it is the cracked edge, go to the next edge
        if (i == crack_edge_IDs(1)) cycle
        ! the two end nodes' coords of edge i
        coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
        coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
        do j = 1, NTESTCRACKPOINT-1
          ! set the temp. 2nd crack point to be the jth test crack point of this edge
          coords_crack(:,2) = (ONE-j*ONE/NTESTCRACKPOINT) * coords_edge(:, 1) + &
          & j * ONE/NTESTCRACKPOINT * coords_edge(:, 2)
          test_line(:)      = coords_crack(:,2) - coords_crack(:,1)
          ! normalize the test_line vector
          call normalize_vect (test_line, is_zero_vect)
          ! if it is too short, then go to the next edge
          if (is_zero_vect) cycle
          ! update projection and 2nd crack_edge_ID & edge crack point
          if (abs(dot_product(test_line,crack_unit_vect)) > projection) then
            projection             = abs(dot_product(test_line,crack_unit_vect))
            crack_edge_IDs(2)      = i
            edge_crack_points(:,2) = coords_crack(:,2)
            n_crack_edge           = 2
          end if
        end do
      end do
      ! if still cannot find the second crack edge, then this elem is wrongly shaped
      if (n_crack_edge == 1) then
        istat = STAT_FAILURE
        emsg  = '2nd crack point cannot be found, element is likely to be &
        & very poorly shaped, crack_elem_centroid2d, global_toolkit_module'
        return
      end if
    ! if none is found, then use two edge midpoints to estimate two crack points
    case (0)
      ! initialize projections
      maxprojection = ZERO
      minprojection = ZERO
      ! start point of test line at centroid
      coords_crack(:,1) = centroid(:)
      do i = 1, nedge
        ! the two end nodes' coords of edge i
        coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
        coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
        do j = 1, NTESTCRACKPOINT-1
          ! set the temp. 2nd crack point to be the jth test crack point of this edge
          coords_crack(:,2) = (ONE-j*ONE/NTESTCRACKPOINT) * coords_edge(:, 1) + &
          & j * ONE/NTESTCRACKPOINT * coords_edge(:, 2)
          test_line(:)      = coords_crack(:,2) - coords_crack(:,1)
          ! normalize the test_line vector
          call normalize_vect (test_line, is_zero_vect)
          ! if it is too short, then go to the next edge
          if (is_zero_vect) cycle
          ! update projections, crack_edge_IDs & edge crack points
          if (dot_product(test_line,crack_unit_vect) > maxprojection) then
            maxprojection          = dot_product(test_line,crack_unit_vect)
            crack_edge_IDs(1)      = i
            edge_crack_points(:,1) = coords_crack(:,2)
          else if (dot_product(test_line,crack_unit_vect) < minprojection) then
            minprojection          = dot_product(test_line,crack_unit_vect)
            crack_edge_IDs(2)      = i
            edge_crack_points(:,2) = coords_crack(:,2)
          end if
        end do
      end do
      ! if still cannot find the second crack edge, then this elem is wrongly shaped
      if (count(crack_edge_IDs>0) /= 2) then
        istat = STAT_FAILURE
        emsg  = '2 crack points cannot be found, element is likely to be &
        & very poorly shaped, crack_elem_centroid2d, global_toolkit_module'
        return
      end if
    ! error
    case default
      istat = STAT_FAILURE
      emsg  = 'unexpected n_crack_edge no. in crack_elem_centroid2d, global_toolkit_module'
      return
    end select


    ! update optional intent out dummy arg only before successful return
    if (present(crack_edge_indices)) crack_edge_indices = crack_edge_IDs

  end subroutine crack_elem_centroid2d



  pure subroutine crack_elem_cracktip2d (cracktip_point, cracktip_edge_index, &
  & nedge, crack_angle, coords, nodes_on_edges, istat, emsg, edge_crack_point,&
  & crack_edge_index)
  ! Purpose:
  ! to find the OTHER cross point between the crack line (starting from existing
  ! crack tip) and another edge of an element, and the index of the edge crossed

  use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                        & ZERO, ONE, HALF, HALFCIRC, PI, SMALLNUM,          &
                        & CROSS_ON_EDGE_ON_CRACK,  CROSS_ON_EDGE_OFF_CRACK, &
                        & NTESTCRACKPOINT

    ! list of dummy args:
    ! - cracktip_point     : coords of the crack tip
    ! - cracktip_edge_index: index of the edge containing the cracktip
    ! - nedge              : no. of edges in this element; = no. of nodes
    ! - crack_angle
    ! - coords             : nodal coordinates of this element
    ! - nodes_on_edges(:,i) : indices of the two end nodes of the edge i
    ! - edge_crack_point   : coords of the crack point on the cracked edge
    ! - crack_edge_index   : index of the cracked edge
    real(DP), intent(in)  :: cracktip_point(2)
    integer,  intent(in)  :: cracktip_edge_index
    integer,  intent(in)  :: nedge
    real(DP), intent(in)  :: crack_angle
    real(DP), intent(in)  :: coords(2, nedge)
    integer,  intent(in)  :: nodes_on_edges(2, nedge)
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    real(DP),                 intent(out) :: edge_crack_point(2)
    integer,        optional, intent(out) :: crack_edge_index

    ! local variables
    ! - theta             : ply angle in radiant
    ! - crack_unit_vector : unit vector along the crack line
    ! - coords_crack      : coordinates of end nodes of crack line segment
    ! - coords_edge       : coordinates of end nodes of an element edge
    ! - cross_point       : coordinates of the cross point btw crack and edge
    ! - n_crack_edge      : no. of cracked edges
    ! - cross_stat        : status of crossing btw a crack and an edge
    ! - crack_edge_ID     : local copy of optional arg crack_edge_index
    real(DP) :: theta, crack_unit_vect(2)
    real(DP) :: coords_crack(2,2), coords_edge(2,2), cross_point(2)
    integer  :: n_crack_edge, cross_stat, crack_edge_ID
    character(len=MSGLENGTH) :: msgloc
    real(DP) :: projection, test_line(2)
    logical  :: is_zero_vect

    integer :: i, j

    ! initialize intent out variables
    istat             = STAT_SUCCESS
    emsg              = ''
    edge_crack_point  = ZERO
    if (present(crack_edge_index)) crack_edge_index = 0
    ! initialize local variables
    theta        = ZERO
    crack_unit_vect = ZERO
    coords_crack = ZERO
    coords_edge  = ZERO
    cross_point  = ZERO
    cross_stat     = 0
    n_crack_edge   = 0
    crack_edge_ID  = 0
    msgloc = ' crack_elem_cracktip2d, global_toolkit_module.'
    projection = ZERO
    test_line  = ZERO
    is_zero_vect = .false.
    i=0; j=0

    ! check validity of inputs: nedge, crack_angle, coords, nodes_on_edges
    ! nedge must be > 0, this check is omitted as if it is not satisfied,
    ! compilor should flag error as it sets the size of coords and nodes_on_edges
    ! crack_angle and coords can take any values
    ! nodes_on_edges must be >= 1, this need to be checked
    if ( cracktip_edge_index < 1 ) then
      istat = STAT_FAILURE
      emsg  = 'cracktip edge index must be >=1,'//trim(msgloc)
      return
    end if
    if ( any(nodes_on_edges < 1) ) then
      istat = STAT_FAILURE
      emsg  = 'end node indices must be >=1,'//trim(msgloc)
      return
    end if

    ! find crack line unit vector
    theta           = crack_angle / HALFCIRC * PI
    crack_unit_vect = [cos(theta), sin(theta)]

    ! set cracktip as node 1 of crack line segment
    coords_crack(:,1) = cracktip_point(:)
    ! set node 2 coordinates of crack line segment
    ! to be a unit distance away from node 1 along the crack line
    coords_crack(:,2) = coords_crack(:,1) + crack_unit_vect(:)



    !**** MAIN CALCULATIONS ****

    ! loop over all edges to find ONE cross point with the crack line
    ! theoretically, this should find EXACTLY ONE cross point, as the
    ! crack line passes the element's cracktip edge and it must cross ANOTHER
    do i = 1, nedge
      ! if edge i is the cracktip edge, cycle to the next one
      if ( i == cracktip_edge_index ) cycle
      ! the two end nodes' coords of edge i
      coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
      coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
      ! zero cross_stat and cross_point for reuse
      cross_stat  = 0
      cross_point = ZERO
      ! check if the edge and crack cross each other
      call edge_crack_cross2d (coords_edge, coords_crack, cross_stat, cross_point)
      ! check cross_stat, update only if cross point is on edge
      if (cross_stat == CROSS_ON_EDGE_ON_CRACK .or. &
      &   cross_stat == CROSS_ON_EDGE_OFF_CRACK) then
          ! the other cracked edge found
          n_crack_edge = 1
          ! store the index of this cracked edge
          crack_edge_ID = i
          ! store the cross point coords of this cracked edge
          edge_crack_point(:) = cross_point(:)
          ! task finished, exit the do loop
          exit
       end if
    end do

    !**** END MAIN CALCULATIONS ****

    !~! check if indeed ONE cracked edge has been found; if not, then the elem
    !~! is likely to be poorly shaped, flag an error, clean up and return
    !~if (n_crack_edge /= 1) then
    !~  istat = STAT_FAILURE
    !~  emsg  = 'Another cracked edge cannot be found, element is likely to be &
    !~  &poorly shaped,'//trim(msgloc)
    !~  ! clean up intent out variable before error exit
    !~  edge_crack_point = ZERO
    !~  return
    !~end if
    
    select case(n_crack_edge)
    ! if the other edge is found, then continue to return
    case (1)
      continue
    ! if the other edge cannot be found, then use an edge midpoint to estimate the other
    case (0)
      ! initialize projection
      projection = ZERO
      ! 1st crack point already found
      coords_crack(:,1) = cracktip_point(:)
      do i = 1, nedge
        ! if it is the cracked edge, go to the next edge
        if (i == cracktip_edge_index) cycle
        ! the two end nodes' coords of edge i
        coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
        coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
        ! loop over all test crack points on this edge
        do j = 1, NTESTCRACKPOINT-1
          ! set the temp. 2nd crack point to be the jth test crack point of this edge
          coords_crack(:,2) = (ONE-j*ONE/NTESTCRACKPOINT) * coords_edge(:, 1) + &
          & j * ONE/NTESTCRACKPOINT * coords_edge(:, 2)
          test_line(:)      = coords_crack(:,2) - coords_crack(:,1)
          ! normalize the test_line vector
          call normalize_vect (test_line, is_zero_vect)
          ! if it is too short, then go to the next edge
          if (is_zero_vect) cycle
          ! update projection and 2nd crack_edge_ID & edge crack point
          if (abs(dot_product(test_line,crack_unit_vect)) > projection) then
            projection          = abs(dot_product(test_line,crack_unit_vect))
            crack_edge_ID       = i
            edge_crack_point(:) = coords_crack(:,2)
            n_crack_edge        = 1
          end if
        end do
      end do
      ! if still cannot find the second crack edge, then this elem is wrongly shaped
      if (n_crack_edge == 0) then
        istat = STAT_FAILURE
        emsg  = 'midpoint of 2nd cracked edge cannot be found, element is likely to be &
        & very poorly shaped,'//trim(msgloc)
        return
      end if
    ! error
    case default
      istat = STAT_FAILURE
      emsg  = 'unexpected n_crack_edge no.,'//trim(msgloc)
      return
    end select

    ! update optional intent out dummy arg only before successful return
    if (present(crack_edge_index)) crack_edge_index = crack_edge_ID

  end subroutine crack_elem_cracktip2d



  pure subroutine partition_quad_elem (NODES_ON_EDGES, crack_edges, &
  & subelem_nodes, istat, emsg, crack_nodes)
  ! Purpose :
  ! partition a quad element into sub elems based on passed in local indices of
  ! its cracked edges
  ! the topology of the quad elem is defined by NODES_ON_EDGES
  ! currently, the algorithm does partition with ONE or TWO cracked edges only
  ! (support for more cracked edges can be included in the future)
  !
  ! **** topology definition of a quad element: ****
  !
  !     4_______10________E3________9________3
  !     |                                    |
  !     |                                    |
  !     |                                    |
  !   11|                                    |8
  !     |                                    |
  !     |                                    |
  !     |                                    |
  !   E4|                                    |E2
  !     |                                    |
  !     |                                    |
  !   12|                                    |7
  !     |                                    |
  !     |                                    |
  !     |____________________________________|
  !     1        5        E1        6        2
  !
  ! **** Each edge is topologically defined by FOUR nodes : ****
  !
  ! nodes on edge E1 ::   <end nodes: 1, 2>   <fl. nodes:  5,  6>
  ! nodes on edge E2 ::   <end nodes: 2, 3>   <fl. nodes:  7,  8>
  ! nodes on edge E3 ::   <end nodes: 3, 4>   <fl. nodes:  9, 10>
  ! nodes on edge E4 ::   <end nodes: 4, 1>   <fl. nodes: 11, 12>
  !
  use parameter_module, only : INT_ALLOC_ARRAY, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

    INTEGER,    PARAMETER :: NNODE = 12, NEDGE = 4

    integer,  intent (in) :: NODES_ON_EDGES(4, NEDGE)
    integer,  intent (in) :: crack_edges(NEDGE)
    type(INT_ALLOC_ARRAY), allocatable, intent(out) :: subelem_nodes(:)
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    integer,        optional, intent(out) :: crack_nodes(4)

    ! local variables
    character(len=MSGLENGTH) :: msgloc
    integer :: n_crackedges, Icrackedge1, Icrackedge2
    integer :: nsub, e1, e2, e3, e4
    integer :: cknodes(4)
    integer :: i, j

    ! initialize intent out and local variables
    istat  = STAT_SUCCESS
    emsg   = ''
    msgloc = ' partition_quad_elem, global_toolkit_module.'
    n_crackedges = 0
    Icrackedge1  = 0
    Icrackedge2  = 0
    nsub = 0
    e1 = 0; e2 = 0; e3 = 0; e4 = 0
    cknodes = 0
    i = 0; j = 0

    ! check input validity
    if ( any (NODES_ON_EDGES < 1) ) then
      istat = STAT_FAILURE
      emsg  = 'incorrect NODES_ON_EDGES in'//trim(msgloc)
      return
    end if
    if ( any(crack_edges < 0) .or. any(crack_edges > NEDGE) ) then
      istat = STAT_FAILURE
      emsg  = 'incorrect crack_edges in'//trim(msgloc)
      return
    end if
    ! check format: lcl ID of cracked edges must be stored first in array elements
    do i = 1, NEDGE-1
        ! when the array elem becomes 0, all the rest must also be 0
        if (crack_edges(i) == 0) then
            ! check the terms after the ith term
            do j = i+1, NEDGE
                ! if any of the following term is NOT 0, flag error and exit program
                if (crack_edges(j) /= 0) then
                  istat = STAT_FAILURE
                  emsg  = "cracked edges' local indices must be stored first in &
                  & crack_edges array in"//trim(msgloc)
                  return
                end if
            end do
        end if
    end do

    ! find no. of cracked edges
    n_crackedges = count (crack_edges > 0)

    select_ncrackedges: select case (n_crackedges)

      ! ***********************************************************************!
      ! one edge cracked, partition this element into a transition element
      ! with three triangular sub elems sharing NODE 3 of the cracked edge
      ! ***********************************************************************!
      case (1) select_ncrackedges

          ! find the index of the cracked edge
          Icrackedge1 = crack_edges(1)

          ! find e1 - e4
          ! e1 - e4: re-index edges to facilitate partitioning domain
          select case(Icrackedge1)
              case (1)
                  e1=1; e2=2; e3=3; e4=4
              case (2)
                  e1=2; e2=3; e3=4; e4=1
              case (3)
                  e1=3; e2=4; e3=1; e4=2
              case (4)
                  e1=4; e2=1; e3=2; e4=3
              case default
                  istat = STAT_FAILURE
                  emsg  = 'wrong 1st broken edge in'//trim(msgloc)
                  return
          end select

          ! nsub: no. of sub elems; here, 3 tri sub elems
          nsub = 3

          ! allocate nsub no. of subelem_nodes
          allocate(subelem_nodes(nsub))
          ! allocate&initialize the internal arrays of these arrays
          do j=1, nsub
            ! allocate connec for 3 nodes of tri sub elems
            allocate(subelem_nodes(j)%array(3))
            ! initialize these arrays
            subelem_nodes(j)%array = 0
          end do

          !:::::::::::::::::::::::::!
          !*** define sub elm 1 ***
          !:::::::::::::::::::::::::!
          ! define its local connec with parent element nodes
          ! top 3 nodes: e1 nodes 3, and e2 nodes 1 & 2
          subelem_nodes(1)%array(1)=NODES_ON_EDGES(3,e1)
          subelem_nodes(1)%array(2)=NODES_ON_EDGES(1,e2)
          subelem_nodes(1)%array(3)=NODES_ON_EDGES(2,e2)
          !:::::::::::::::::::::::::!
          !*** define sub elm 2 ***
          !:::::::::::::::::::::::::!
          ! define its local connec with parent element nodes
          ! top 3 nodes: e2 node 3, and e3 nodes 1 & 2
          subelem_nodes(2)%array(1)=NODES_ON_EDGES(3,e1)
          subelem_nodes(2)%array(2)=NODES_ON_EDGES(1,e3)
          subelem_nodes(2)%array(3)=NODES_ON_EDGES(2,e3)
          !:::::::::::::::::::::::::!
          !*** define sub elm 3 ***
          !:::::::::::::::::::::::::!
          ! define its local connec with parent element nodes
          ! top 3 nodes: e4 nodes 1 & 2, and e1 node 3
          subelem_nodes(3)%array(1)=NODES_ON_EDGES(1,e4)
          subelem_nodes(3)%array(2)=NODES_ON_EDGES(2,e4)
          subelem_nodes(3)%array(3)=NODES_ON_EDGES(3,e1)

          ! ** NOTE : ONLY NODE 3 of the cracked edge is used for this partition


      ! ***********************************************************************!
      ! TWO edge cracked, partition this element into two quad sub elems or
      ! four triangular sub elems
      ! ***********************************************************************!
      case (2) select_ncrackedges
      !- two edges cracked

          ! find the indices of the two broken edges, sorted in ascending order
          Icrackedge1 = min( crack_edges(1), crack_edges(2) )
          Icrackedge2 = max( crack_edges(1), crack_edges(2) )

          ! Icrackedge1 must be between 1 to 3, and Icrackedge2 between 2 to 4,
          ! with Icrackedge2 > Icrackedge1
          if (Icrackedge1<1 .or. Icrackedge1>3 .or. &
          &   Icrackedge2<2 .or. Icrackedge2>4 .or. &
          &   Icrackedge2 <= Icrackedge1) then
            istat = STAT_FAILURE
            emsg = 'wrong crack edge indices, &
            & case n_crackedges=2, in'//trim(msgloc)
            return
          end if

          ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
          ! determine partition based on the indices of the two broken edges
          ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
          ! find nsub and e1 - e4
          ! nsub: no. of sub domains
          ! e1 - e4: re-index edges to facilitate partitioning domain
          select case(Icrackedge1)
            case(1)
            ! 1st cracked edge is lcl edge 1, decide on nsub and e1-e4 based on
            ! the lcl ID of the 2nd cracked edge
                select case (Icrackedge2)
                  case(2)
                    nsub=4
                    e1=1; e2=2; e3=3; e4=4
                  case(3)
                    nsub=2
                    e1=1; e2=2; e3=3; e4=4
                  case(4)
                    nsub=4
                    e1=4; e2=1; e3=2; e4=3
                  case default
                    istat = STAT_FAILURE
                    emsg = 'wrong 2nd broken edge in'//trim(msgloc)
                    return
                end select
            case(2)
                select case (Icrackedge2)
                  case(3)
                    nsub=4
                    e1=2; e2=3; e3=4; e4=1
                  case(4)
                    nsub=2
                    e1=2; e2=3; e3=4; e4=1
                  case default
                    istat = STAT_FAILURE
                    emsg = 'wrong 2nd broken edge in'//trim(msgloc)
                    return
                end select
            case(3)
                if(Icrackedge2==4) then
                    nsub=4
                    e1=3; e2=4; e3=1; e4=2
                else
                    istat = STAT_FAILURE
                    emsg = 'wrong 2nd broken edge in'//trim(msgloc)
                    return
                end if
            case default
                istat = STAT_FAILURE
                emsg = 'wrong 1st broken edge in'//trim(msgloc)
                return
          end select

          ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
          ! form sub elements based on nsub values
          ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
          select_nsub: select case(nsub)

            ! :::::::::::::::::::::::::::::::::::::::::!
            ! two quad subdomains
            ! :::::::::::::::::::::::::::::::::::::::::!
            case(2) select_nsub
                !:::::::::::::::::::::::::!
                !*** prepare arrays ***
                !:::::::::::::::::::::::::!
                ! allocate nsub no. of subelem_nodes
                allocate(subelem_nodes(nsub))
                ! allocate&initialize the internal arrays of these arrays
                do j=1, nsub
                  ! allocate connec for the 4 nodes of quad sub elems
                  allocate(subelem_nodes(j)%array(4))
                  ! initialize these arrays
                  subelem_nodes(j)%array = 0
                end do
                !:::::::::::::::::::::::::!
                !*** define sub elm 1 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 4 nodes: e1 nodes 1 & 3, and e3 nodes 4 & 2
                subelem_nodes(1)%array(1)=NODES_ON_EDGES(1,e1)
                subelem_nodes(1)%array(2)=NODES_ON_EDGES(3,e1)
                subelem_nodes(1)%array(3)=NODES_ON_EDGES(4,e3)
                subelem_nodes(1)%array(4)=NODES_ON_EDGES(2,e3)
                !:::::::::::::::::::::::::!
                !*** define sub elm 2 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 4 nodes: e1 nodes 4 & 2, and e3 nodes 1 & 3
                subelem_nodes(2)%array(1)=NODES_ON_EDGES(4,e1)
                subelem_nodes(2)%array(2)=NODES_ON_EDGES(2,e1)
                subelem_nodes(2)%array(3)=NODES_ON_EDGES(1,e3)
                subelem_nodes(2)%array(4)=NODES_ON_EDGES(3,e3)
                !:::::::::::::::::::::::::!
                !*** define cknodes ***
                !:::::::::::::::::::::::::!
                cknodes(1) = NODES_ON_EDGES(4,e1)
                cknodes(2) = NODES_ON_EDGES(3,e3)
                cknodes(3) = NODES_ON_EDGES(4,e3)
                cknodes(4) = NODES_ON_EDGES(3,e1)

            ! :::::::::::::::::::::::::::::::::::::::::!
            ! four tri subdomains
            ! :::::::::::::::::::::::::::::::::::::::::!
            case(4) select_nsub
                !:::::::::::::::::::::::::!
                !*** prepare arrays ***
                !:::::::::::::::::::::::::!
                ! allocate nsub no. of subelem_nodes
                allocate(subelem_nodes(nsub))
                ! allocate&initialize the internal arrays of these arrays
                do j=1, nsub
                  ! allocate connec for 3 nodes of tri sub elems
                  allocate(subelem_nodes(j)%array(3))
                  ! initialize these arrays
                  subelem_nodes(j)%array = 0
                end do
                !:::::::::::::::::::::::::!
                !*** define sub elm 1 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 3 nodes: e1 nodes 4, and e2 nodes 1 & 3
                subelem_nodes(1)%array(1)=NODES_ON_EDGES(4,e1)
                subelem_nodes(1)%array(2)=NODES_ON_EDGES(1,e2)
                subelem_nodes(1)%array(3)=NODES_ON_EDGES(3,e2)
                !:::::::::::::::::::::::::!
                !*** define sub elm 2 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 3 nodes: e2 node 4, and e3 nodes 1 & 2
                subelem_nodes(2)%array(1)=NODES_ON_EDGES(4,e2)
                subelem_nodes(2)%array(2)=NODES_ON_EDGES(1,e3)
                subelem_nodes(2)%array(3)=NODES_ON_EDGES(2,e3)
                !:::::::::::::::::::::::::!
                !*** define sub elm 3 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 3 nodes: e4 nodes 1 & 2, and e1 node 3
                subelem_nodes(3)%array(1)=NODES_ON_EDGES(1,e4)
                subelem_nodes(3)%array(2)=NODES_ON_EDGES(2,e4)
                subelem_nodes(3)%array(3)=NODES_ON_EDGES(3,e1)
                !:::::::::::::::::::::::::!
                !*** define sub elm 4 ***
                !:::::::::::::::::::::::::!
                ! define its local connec with parent element nodes
                ! top 3 nodes: e1 node 3, e2 node 4 and e3 node 2
                subelem_nodes(4)%array(1)=NODES_ON_EDGES(3,e1)
                subelem_nodes(4)%array(2)=NODES_ON_EDGES(4,e2)
                subelem_nodes(4)%array(3)=NODES_ON_EDGES(2,e3)
                !:::::::::::::::::::::::::!
                !*** define cknodes ***
                !:::::::::::::::::::::::::!
                cknodes(1) = NODES_ON_EDGES(4,e1)
                cknodes(2) = NODES_ON_EDGES(3,e2)
                cknodes(3) = NODES_ON_EDGES(4,e2)
                cknodes(4) = NODES_ON_EDGES(3,e1)
                
            ! :::::::::::::::::::::::::::::::::::::::::!
            ! unsupported no. of sub elems, ERROR
            ! :::::::::::::::::::::::::::::::::::::::::!
            case default select_nsub
                istat = STAT_FAILURE
                emsg = 'wrong nsub in'//trim(msgloc)
                return
          end select select_nsub


      ! ***********************************************************************!
      ! Unsupported number of cracked edge
      ! return with error
      ! ***********************************************************************!
      case default select_ncrackedges
          istat = STAT_FAILURE
          emsg = 'unexpected n_crackedges value in'//trim(msgloc)
          return

      end select select_ncrackedges
      
      
      if (present(crack_nodes)) crack_nodes = cknodes


  end subroutine partition_quad_elem



  pure subroutine assembleKF (Kmat, Fvec, Ki, Fi, cnc, NDIM, istat, emsg)
  ! Purpose:
  ! to assemble elem/subelem K and F to system/xelem K and F
  ! variable sizes of dummy args have to be allowed here.
  ! explicit checking of dummy args' sizes must be performed
  ! with informative error messages flagged if unexpected sizes
  ! are encountered
  !
  ! the inputs must safisty the following conditions:
  ! - Kmat is square matrix with size = size of Fvec
  ! - Ki   is square matrix with size = size of Fi
  ! - size of cnc * NDIM = size of Fi
  ! - min element of cnc must > 0
  ! - max element of cnc * NDIM must < size of Fvec

  use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

    real(DP), intent(inout) :: Kmat(:,:), Fvec(:)
    real(DP), intent(in)    :: Ki(:,:),   Fi(:)
    integer,  intent(in)    :: cnc(:)
    integer,  intent(in)    :: NDIM
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg

    ! local variable
    integer,  allocatable :: dofcnc(:)
    integer :: i, j, l, nF, nFi, ncnc, mincnc, maxcnc

    ! initialize intent out and local variables
    istat  = STAT_SUCCESS
    emsg   = ''
    i      = 0
    j      = 0
    l      = 0
    nF     = 0
    nFi    = 0
    ncnc   = 0
    mincnc = 0
    maxcnc = 0

    nF     = size(Fvec)
    nFi    = size(Fi)
    ncnc   = size(cnc)
    mincnc = minval(cnc)
    maxcnc = maxval(cnc)

    ! check validity of inputs
    if ( any (shape(Kmat) /= [nF,nF]) ) then
      istat = STAT_FAILURE
      emsg  = 'Kmat shape is incorrect, assembleKF, global_toolkit_module'
    else if ( any (shape(Ki) /= [nFi,nFi]) ) then
      istat = STAT_FAILURE
      emsg  = 'Ki shape is incorrect, assembleKF, global_toolkit_module'
    else if ( ncnc * NDIM /= nFi ) then
      istat = STAT_FAILURE
      emsg  = 'cnc size is incorrect, assembleKF, global_toolkit_module'
    else if ( mincnc < 1 ) then
      istat = STAT_FAILURE
      emsg  = 'cnc min element <1, assembleKF, global_toolkit_module'
    else if ( maxcnc * NDIM > nF ) then
      istat = STAT_FAILURE
      emsg  = 'cnc max element is too large for Kmat, assembleKF, &
      &global_toolkit_module'
    end if

    if (istat == STAT_FAILURE) return

    ! proceed with the assembly only when all dummy arg sizes are consistent
    
    ! prepare the d.o.f connectivity array, dofcnc
    allocate(dofcnc(nFi)); dofcnc = 0
    
    ! loop over no. of nodes in this elem to form its dofcnc
    do j = 1, ncnc
      do l = 1, NDIM
        ! dof indices of the jth node of this elem
        dofcnc( (j-1) * NDIM + l ) = ( cnc(j) - 1 ) * NDIM + l
      end do
    end do
    
    do i = 1, nFi
      do j = 1, nFi
          Kmat(dofcnc(j),dofcnc(i)) = Kmat(dofcnc(j),dofcnc(i)) + Ki(j,i)
      end do
      Fvec(dofcnc(i)) = Fvec(dofcnc(i)) + Fi(i)
    end do

  end subroutine assembleKF





end module global_toolkit_module