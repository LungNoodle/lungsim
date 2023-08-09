module capillaryflow
!*Brief Description:* This module handles all microcirculatory blood flow.
!
!*LICENSE:*
!
!
!*Contributor(s):* Kelly Burrowes, Alys Clark
!
!*Full Description:*
!
  !This module handles all microcirculatory blood flow.

  use arrays
  use diagnostics
  use indices
  use math_utilities
  use other_consts
  use precision
  use solve

  implicit none


  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public cap_flow_ladder,cap_flow_admit

contains
!
!###################################################################################
!
  subroutine cap_flow_ladder(ne,LPM_FUNC,Lin,Lout,Pin,Pout,Ppl,&
    Q01,R_in,R_out,x,y,z,OUTPUT_PERFUSION)
!*cap_flow_ladder:* Uses the ladder model to solve for perfusion in the acinus.
! It uses a symmetric, bifuracting arteriole and venule tree with N generations
! and calls a function to calculate resistance across a capillary sheet at each
! generation (CAP_FLOW_SHEET.f).

!###  -----------------------INPUT-------------------------
!###  The input to this subroutine is:
!###  ne=element number
!###  Pin= Pressure into the acinus
!###  Pout=Pressure out of the acinus
!###  Ppl= Pleural pressure

!###  ----------------------OUTPUT----------------------------
!###  Important output to large vessel models
!###  LPM_FUNC= Resistance across the acinus:

!###  This needs to be fed back into the large vessel model. Then
!###  Pin-Pout=LPM_FUNC*Q is solved for Pin, Pout and Q as part
!###  that system. Pin and Pout can then be fed back iteratively
!###  into this subroutine...
!###
!###  In addition this subroutine outputs:
!###  Pressure at each each vessel interection and flow, resistance
!###  and RBC transit times through each capillary element.

!###  UNITS. The units that are used here are m.

    integer, intent(in) :: ne
    real(dp), intent(inout) :: LPM_FUNC
    real(dp):: Pin,Pout,Q01,R_in,R_out,x,y,z,Lin,Lout,Ppl
    logical, intent(in) :: OUTPUT_PERFUSION

    type(capillary_bf_parameters) :: cap_param

    !    Local variables
    integer :: MatrixSize,NonZeros,submatrixsize,ngen,i
    real(dp) :: area,Q01_mthrees,sheet_number
    character(len=60) :: sub_name

    sub_name = 'cap_flow_ladder'
    call enter_exit(sub_name,1)

!     Number of non-zero entries in solution matrix.
      NonZeros=3
      do i=2,cap_param%num_symm_gen
         NonZeros=NonZeros+4*i+10
      enddo
!     The size of the solution matrix (number of unknown pressures and flows)
      MatrixSize=5*cap_param%num_symm_gen-3
!!     The number of unknown pressures
      submatrixsize=4*cap_param%num_symm_gen-4
      ngen=cap_param%num_symm_gen

!...  ---INITIALISATION
!...  The input Q01 gives us an estimate for flow into the acinus from the large
!...  vessel model.
!...  This is in mm^3/s and needs to be converted to m^3/s to use in calculating
!...  arteriole and venule resistance
      Q01_mthrees=Q01/1.d9 !mm3/s->m3/s
!     Sheet area (unscaled):
!...  We define a sheet area for input into the capillary model.
!...  This area is at full inflation and will be scaled within CAP_FLOW_SHEET
!...  Area of an individual sheet
      sheet_number=0
      DO i=1,cap_param%num_symm_gen
         sheet_number=sheet_number+2.d0**i
      ENDDO
      area=cap_param%total_cap_area/(sheet_number*num_units) !m^2
!
!!  ---CALL THE FUNCTIONS THAT CALCULATE THE FLOW ACROSS THE LADDER FOR A GIVEN PRESSURE DROP--
      call evaluate_ladder(ne,NonZeros,MatrixSize,submatrixsize,ngen,&
      area,Lin,Lout,Pin,Pout,Ppl,R_in,R_out,Q01_mthrees,x,y,z,&
      OUTPUT_PERFUSION)

!     ---FINAL FUNCTION OUTPUT (Resistance across ladder)---
!...  This takes difference between the inlet and outlet pressures
!...  (Pin and Pout) and divides by an updated flow (Q01_mthrees)
!...  to give updated resistance across the ladder. This feeds back
!...  to the large vessel model.
      LPM_FUNC=(Pin-Pout)/(Q01_mthrees*1000.d0**3) !Pa.s/m^3->pa.s/mm^3
    call enter_exit(sub_name,2)

  end subroutine cap_flow_ladder
!
!######################################################################################################
!
!*evaluate_ladder:* Sets up and solves matrix equations for ladder model.
  subroutine evaluate_ladder(ne,NonZeros,MatrixSize,submatrixsize,ngen,&
      area,L_in,L_out,Pin,Pout,Ppl,R_in,R_out,Q01_mthrees,x,y,z, &
      OUTPUT_PERFUSION)

    type(capillary_bf_parameters) :: cap_param

    integer :: ne,NonZeros,MatrixSize,submatrixsize,ngen
    real(dp) :: area,Pin,Pout,Ppl,Q01_mthrees,x,y,z,R_in,R_out,L_in,L_out
    logical :: OUTPUT_PERFUSION
    real(dp) :: area_scale,length_scale,alpha_c

    ! Local variables
    integer :: i,iter,j,gen,zone,num_sheet
    integer, allocatable :: SparseCol(:)
    integer, allocatable ::SparseRow(:)
    real(dp) :: area_new,ErrorEstimate,Hart,Hven,Pin_SHEET,Pout_SHEET
    real(dp), allocatable :: Pressure(:)
    real(dp) ::  Q_c,Qtot,Qgen
    real(dp),allocatable :: Q_sheet(:)
    real(dp),allocatable :: RHS(:)
    real(dp) :: RBC_TT,Rtot,SHEET_RES
    real(dp),allocatable :: Solution(:)
    real(dp),allocatable :: SolutionLast(:)
    real(dp),allocatable :: SparseVal(:)
    real(dp) :: TOTAL_CAP_VOL,TOTAL_SHEET_H,TOTAL_SHEET_SA,&
         recruited,TT_TOTAL,R_upstream,R_downstream
    real(dp),allocatable :: l_a(:),rad_a(:),l_v(:),rad_v(:),mu_app(:)
     integer SOLVER_FLAG
    character(len=60) :: sub_name
    integer :: AllocateStatus

    sub_name = 'evaluate_ladder'
    call enter_exit(sub_name,1)

    allocate (Pressure(submatrixsize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SparseCol(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SparseVal(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (Solution(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SolutionLast(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (l_a(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (rad_a(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (l_v(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (rad_v(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (Q_Sheet(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
     allocate (mu_app(ngen), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"

    Pressure=0.0_dp
    SparseCol=0
    SparseRow=0
    SparseVal=0.0_dp
    Solution=0.0_dp
    SolutionLast=0.0_dp

!!...  Initial guess for pressure distribution lets say all arterial pressures are the same
!!...  and all the venous pressures are the same solution appears independent of this.
      DO i=1,cap_param%num_symm_gen-1
          Pressure(4*i-3)=1000.0_dp ! Pa
          Pressure(4*i-2)=1000.0_dp
          Pressure(4*i-1)=100.0_dp
          Pressure(4*i)=100.0_dp
      ENDDO


!###  INITIAL SOLUTION GUESS
      DO i=1,submatrixsize
        solution(i)=Pressure(i)
      ENDDO
      DO i=submatrixsize+1,matrixSize-1
       solution(i)=Q01_mthrees/2**cap_param%num_symm_gen
      ENDDO
      Solution(Matrixsize)=Q01_mthrees
!###  INITIALISE SOLUTIONLAST
      DO j=1,MatrixSize
        SolutionLast(j)=Solution(j)
      ENDDO

!### INPUT TO THE LADDER MODEL THAT IS INDEPENDENT OF ITERATION
      CALL LADDERSOL_MATRIX(NonZeros,MatrixSize,submatrixsize,&
         SparseCol,SparseRow,SparseVal,RHS,Pin,Pout)

      CALL cap_specific_parameters(ne,Ppl,alpha_c,area_scale,length_scale,l_a,rad_a,l_v,rad_v,ngen,&
        mu_app,R_in,R_out,L_in,L_out)
!### ITERATIVE LOOP
      iter=0
      ErrorEstimate=1.d10
      DO WHILE(ErrorEstimate.GT.1.0d-9.AND.iter.LT.100)
        iter=iter+1
!...  CALCULATE RESISTANCE GIVEN CURRENT PRESSURE AND FLOW - THEN UPDATE
!...  SparseVal- THese are the only elements of the solution matrix that need
!.... iteratively updating

        CALL POPULATE_MATRIX_LADDER(ne,NonZeros,submatrixsize,ngen,area,alpha_c,&
          area_scale,length_scale,mu_app,Ppl,Pin,Pout,Pressure,&
         Q01_mthrees,Q_sheet,SparseVal,l_a,rad_a,l_v,rad_v)


      call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
         Solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)

!        call BICGSTAB_LinSolv(MatrixSize,NonZeros,RHS,Solution,SparseCol,&
!             SparseRow,SparseVal,1.d-9,1000) !NB/ 1.d-9 = convergence tol, 1000 = max solver iterations



         DO j=1,submatrixsize
            Pressure(j)=Solution(j)
         ENDDO
         Q01_mthrees=Solution(MatrixSize)
!     Estimating Error in solution
      ErrorEstimate=0.d0
        DO i=1,MatrixSize
          ErrorEstimate=ErrorEstimate+&
          DABS((Solution(i)-SolutionLast(i))**2.d0&
          /Solution(i)**2.d0)
          SolutionLast(i)=Solution(i)
        ENDDO
      ErrorEstimate=ErrorEstimate/MatrixSize

      ENDDO

      Qtot=0
      Do i=1,cap_param%num_symm_gen
        Qtot=Qtot+Q_sheet(i)*2.d0**i
      ENDDO
      Rtot=(Pin-Pout)/Q01_mthrees
!
       IF(OUTPUT_PERFUSION)THEN
!###  GET SOLUTIONS TO WRITE TO FILE
        TOTAL_CAP_VOL=0.d0
        TOTAL_SHEET_SA=0.d0
        TT_TOTAL=0.d0
        TOTAL_SHEET_H=0.d0
        num_sheet=0
!... Ladder output
        DO i=1,cap_param%num_symm_gen-1
           gen=i
           Pin_sheet=Pressure(4*i-3)
           Pout_sheet=Pressure(4*i-1)
!...     calulate resistance and transit time through a single capillary
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,area,area_new,alpha_c,area_scale,&
        length_scale,recruited)
           num_sheet=num_sheet+2**gen
           Qgen=Q_c*2.d0**i
         WRITE(10,&
              '(I6,X,3(F9.2,X),I6,X,4(F8.2,X),4(F8.5,X),&
         2(F10.2,X),3(F8.4,X),I6,X,2(F10.5,X),2(F8.4,X),&
         (F10.2,X))')&
         ne,x,y,z,gen,Pin,Pin_sheet,Pout_sheet,Pout,Qtot*1.d9,&
         Qgen*1.d9,Q_c*1.d9,(Hart-Hven)*1.d6,SHEET_RES/1000.d0**3.d0,&
         Rtot/1000.d0**3.d0,RBC_tt,Hart*1.d6,Hven*1.d6,zone,&
         (Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen),&
         area_new*1.d6*(2.d0**gen),recruited

           TOTAL_CAP_VOL=TOTAL_CAP_VOL+(Hart/2.d0+Hven/2.d0)*&
                area_new*1.d9*(2.d0**gen)
           TOTAL_SHEET_H=TOTAL_SHEET_H+(Hart/2.d0+Hven/2.d0)*(2.d0**gen)*1.d6
          TOTAL_SHEET_SA=TOTAL_SHEET_SA+area_new*1.d6*(2.d0**gen)
          TT_TOTAL=TT_TOTAL+RBC_tt*(2.d0**gen)
        ENDDO

        gen=cap_param%num_symm_gen
        Pin_sheet=Pressure(4*cap_param%num_symm_gen-6) !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*cap_param%num_symm_gen-4) !%pressure out of final capillary sheets
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,area,area_new,alpha_c,area_scale,&
        length_scale,recruited)
           Qgen=Q_C*2.d0**cap_param%num_symm_gen
         WRITE(10,&
        '(I6,X,3(F9.2,X),I6,X,4(F8.2,X),4(F8.5,X),&
        2(F10.2,X),3(F8.4,X),I6,X,2(F10.5,X),3(F8.4,X))')&
         ne,x,y,z,gen,Pin,Pin_sheet,Pout_sheet,Pout,Qtot*1.d9,&
         Qgen*1.d9,Q_c*1.d9,(Hart-Hven)*1.d6,SHEET_RES/1000.d0**3.d0,&
         Rtot/1000.d0**3.d0,RBC_tt,Hart*1.d6,Hven*1.d6,zone,&
         R_upstream,R_downstream,&
         (Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen),&
         area_new*1.d6*(2.d0**gen),recruited

          TOTAL_CAP_VOL=TOTAL_CAP_VOL&
            +(Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen)
           TOTAL_SHEET_H=TOTAL_SHEET_H&
            +(Hart/2.d0+Hven/2.d0)*(2.d0**gen)*1.d6
          TOTAL_SHEET_H=TOTAL_SHEET_H/num_sheet
          TOTAL_SHEET_SA=TOTAL_SHEET_SA&
            +area_new*1.d6*(2.d0**gen)
            TT_TOTAL=TT_TOTAL+RBC_tt*(2.d0**gen)
            TT_TOTAL=TT_TOTAL/num_sheet

!... General output
! ne=1  |  x=2  |  y=3  |  z=4  | Pin=5 Pa |Pout=6 Pa | Qtot=7 mm^3/s |sum Qsheet=8 mm^3 /s|
! Rtot=9 Pa/mm^3 | Blood_vol=10 mm^3| sheet_area= 11 mm^2 | ave_TT=12 s |ave_H=13 um |Ppl=14 Pa
          WRITE(20,&
        '(I6,X,5(F9.2,X),2(F8.5,X),F10.2,X,F8.4,X,F10.4,X,F10.3,X,F8.4,X,F9.4,X)') &
         ne,x,y,z,Pin,Pout,Q01_mthrees*1.d9,Qtot*1.d9,Rtot/1000.d0**3.d0,&
         TOTAL_CAP_VOL,TOTAL_SHEET_SA,TT_TOTAL,TOTAL_SHEET_H,Ppl
        ENDIF
!
    deallocate (SparseCol, STAT = AllocateStatus)
    deallocate (SparseRow, STAT = AllocateStatus)
    deallocate (SparseRow, STAT = AllocateStatus)
    deallocate (Solution, STAT = AllocateStatus)
    deallocate (SolutionLast, STAT = AllocateStatus)
    deallocate (Pressure, STAT = AllocateStatus)
    deallocate (RHS, STAT = AllocateStatus)
    deallocate (l_a, STAT = AllocateStatus)
    deallocate (l_v, STAT = AllocateStatus)
    deallocate (rad_a, STAT = AllocateStatus)
    deallocate (rad_v, STAT = AllocateStatus)
    deallocate (Q_sheet, STAT = AllocateStatus)
    deallocate (mu_app, STAT = AllocateStatus)


    call enter_exit(sub_name,2)

  end subroutine evaluate_ladder
!
!######################################################################################################
!
subroutine laddersol_matrix(NonZeros,MatrixSize,submatrixsize,&
      SparseCol,SparseRow,SparseVal,RHS,Pin,Pout)
!*laddersol_matrix:*Sets up ladder matrix entries that are independent of iteration.

    type(capillary_bf_parameters) :: cap_param

    integer :: NonZeros,MatrixSize,submatrixsize,SparseRow(MatrixSize+1),&
       SparseCol(Nonzeros)
    real(dp) :: SparseVal(Nonzeros),RHS(MatrixSize),Pin,Pout
    !    Local variables
    integer :: count1,count,i,j
    character(len=60) :: sub_name

    sub_name = 'laddesol_martrix'
    call enter_exit(sub_name,1)

!  Define the vector Ap (for UMFPACK) this vector starts at zero and then each entry is the cumulative total of non-zero entries as you step through columns from L to right
      SparseRow(1)=1
      SparseRow(2)=3
      SparseRow(3)=7
      SparseRow(4)=9
      SparseRow(5)=13
      DO i=2,cap_param%num_symm_gen-1
         SparseRow(6+4*(i-2))=SparseRow(6+4*(i-2)-1)+2+i
         SparseRow(7+4*(i-2))=SparseRow(6+4*(i-2))+3+i
         SparseRow(8+4*(i-2))=SparseRow(7+4*(i-2))+2+i
         SparseRow(9+4*(i-2))=SparseRow(8+4*(i-2))+3+i
      ENDDO
      DO i=1,cap_param%num_symm_gen
        SparseRow(4*cap_param%num_symm_gen-3+i)=SparseRow(submatrixsize+i)+3
      ENDDO
      SparseRow(MatrixSize+1)=SparseRow(MatrixSize)+cap_param%num_symm_gen+1

!  DEFINE THE RHS VECTOR (size=Matrix size)- These are the BCS
      DO count1=1,MatrixSize
         RHS(count1)=0.d0
      ENDDO
      RHS(1)=-Pin
      RHS(3)=Pout

! Define the column indices and Ax values for UMFPACK. These step through the rowss and give the column index and the value. (start at column 1)

      !...  ENTRIES THAT ARE INDEPENDENT OF ITERATION
      DO i=1,Nonzeros
         SparseCol(i)=0
         SparseVal(i)=0.d0
      ENDDO
      !...  CONSERVATION OF FLOW
      SparseCol(NonZeros)=MatrixSize
      SparseVal(NonZeros)=1.d0
      DO i=1,cap_param%num_symm_gen
         SparseCol(NonZeros-i)=MatrixSize-i
         SparseVal(NonZeros-i)=-2.d0**(cap_param%num_symm_gen+1-i)
         SparseCol(NonZeros-cap_param%num_symm_gen-1-3*(i-1))=MatrixSize-i
      ENDDO
      !...  CAPILLARY CONNECTIONS
      !...  Prior to final generation..
      DO i=1,cap_param%num_symm_gen-1
        SparseCol(NonZeros-4*cap_param%num_symm_gen+3*(i-1))=1+4*(i-1)
        SparseVal(NonZeros-4*cap_param%num_symm_gen+3*(i-1))=1.d0
        SparseCol(NonZeros-4*cap_param%num_symm_gen+3*(i-1)+1)=3+4*(i-1)
        SparseVal(NonZeros-4*cap_param%num_symm_gen+3*(i-1)+1)=-1.d0
      ENDDO
      !...  FINAL GENERATION
      SparseCol(NonZeros-cap_param%num_symm_gen-2)=MatrixSize-cap_param%num_symm_gen-1
      SparseVal(NonZeros-cap_param%num_symm_gen-2)=-1.d0
      SparseCol(NonZeros-cap_param%num_symm_gen-3)=MatrixSize-cap_param%num_symm_gen-3
      SparseVal(NonZeros-cap_param%num_symm_gen-3)=1.d0
      !...  First generation
      SparseCol(1)=1
      Sparseval(1)=-1.d0
      SparseCol(2)=MatrixSize
      SparseCol(3)=1
      SparseVal(3)=1.d0
      SparseCol(4)=2
      SparseVal(4)=-1.d0
      SparseCol(5)=submatrixsize+1
      SparseCol(6)=MatrixSize
      SparseCol(7)=3
      SparseVal(7)=1.d0
      SparseCol(8)=MatrixSize
      SparseCol(9)=3
      SparseVal(9)=-1.d0
      SparseCol(10)=4
      SparseVal(10)=1.d0
      SparseCol(11)=submatrixsize+1
      SparseCol(12)=MatrixSize
      count=12
      DO i=2,cap_param%num_symm_gen-1
      SparseCol(count+1)=4*i-6
      SparseVal(count+1)=1.d0
      SparseCol(count+2)=4*i-3
      SparseVal(count+2)=-1.d0
        DO j=2,i
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+2
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=5+4*(i-2)
      SparseVal(count+1)=1.d0
      SparseCol(count+2)=6+4*(i-2)
      SparseVal(count+2)=-1.d0
        DO j=2,i+1
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+3
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=4*i-4
      SparseVal(count+1)=-1.d0
      SparseCol(count+2)=4*i-1
      SparseVal(count+2)=1.d0
         DO j=2,i
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+2
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=4*i-1
      SparseVal(count+1)=-1.d0
      SparseCol(count+2)=4*i
      SparseVal(count+2)=1.d0
        DO j=2,i+1
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+3
      SparseCol(count)=MatrixSize
      ENDDO

    call enter_exit(sub_name,2)

end subroutine laddersol_matrix
!
!######################################################################################################
!
subroutine populate_matrix_ladder(ne,NonZeros,submatrixsize,ngen,area,alpha_c,&
       area_scale,length_scale,mu_app,Ppl,Pin,Pout,Pressure,&
      Q01_mthrees,Q_sheet,SparseVal,l_a,rad_a,l_v,rad_v)
!*populate_matrix_ladder:*Sets up ladder matrix entries that are NOT independent of iteration.

    type(capillary_bf_parameters) :: cap_param

    integer :: ne,NonZeros,submatrixsize,ngen
    real(dp) :: area,mu_app(ngen),Pin,Pout,Ppl,Pressure(submatrixsize),Q01_mthrees,&
            Q_sheet(ngen),SparseVal(NonZeros),alpha_c,&
            area_scale,length_scale,l_a(ngen),rad_a(ngen),l_v(ngen),rad_v(ngen)



    !    Local variables
    integer :: gen,count,j,zone
    real(dp) :: Q,P_exta,P_extv,radupdate,R_art1,R_art2,&
        R_ven1,R_ven2,SHEET_RES,Q_c,R_sheet(ngen),Q_gen,Hart,&
        Hven,RBC_TT,Pin_sheet,Pout_sheet,area_new,test,&
        recruited,volume_vessels

    character(len=60) :: sub_name

    sub_name = 'populate_matrix_ladder'
    call enter_exit(sub_name,1)
!   RESISTANCE CALCULATIONS STEPPING THROUGH THE GENERATIONS
!   AND CONSERVATION OF FLOW!!
!...  initialising vessel volumes
      volume_vessels=0.d0

!...  Previous iterations estimate for total flow through the system
      !ALYS: at the moment we aren't iterating so use Q01
       Q=Q01_mthrees
       radupdate=0.d0
      DO gen=1,cap_param%num_symm_gen-1
!...    FIRST HALF OF ARTERIOLE
!...    Update radius of arteriole based on inlet pressure
         IF(rad_a(gen).LT.100.d-6) THEN
           P_exta=cap_param%Palv ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_exta=-Ppl
         ENDIF
         IF ((Pressure(4*gen-3)-P_exta).LE.cap_param%Pub_a_v)THEN
!            radupdate=rad_a(gen)*(1+(Pressure(4*gen-3)-P_exta)*cap_param%alpha_a)
            radupdate=rad_a(gen)+cap_param%alpha_a*(Pressure(4*gen-3)-P_exta)*&
            (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-3)-P_exta)/(2.d0-2.d0*cap_param%num_symm_gen)
         ELSE
!            radupdate=rad_a(gen)*(1+cap_param%Pub_a_v*cap_param%alpha_a)
            radupdate=rad_a(gen)+cap_param%alpha_a*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
            /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
            /(2.d0-2.d0*cap_param%num_symm_gen)
         ENDIF
!         IF(nj_hypoxia.NE.0) THEN
!           IF(rad_a(gen).LE.0.25d0)THEN
!              radupdate=radupdate*k_factor
!           ENDIF
!         ENDIF
         volume_vessels=volume_vessels+&
          (2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0)

!...  Calculate Poiseuille resistance in first half of arteriole - (only
!...  half total generation length)
           R_art1=(8.d0*mu_app(gen)*L_a(gen)/2.d0)/(pi*radupdate**4.d0)

!...    FIRST HALF OF VENULE
!...    Update radius of venule based on inlet pressure NB: May need to
!...     update this later to give an average of inlet and outlet radii
           IF (Pressure(4*gen-1)-cap_param%Palv.LE.cap_param%Pub_a_v)THEN
           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-cap_param%Palv)*cap_param%alpha_v)
           ELSE
              radupdate=rad_v(gen)*(1+cap_param%Pub_a_v*cap_param%alpha_v)
           ENDIF
         IF(rad_v(gen).LT.100.d-6) THEN
           P_extv=cap_param%Palv ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_extv=-Ppl
         ENDIF
         IF((Pressure(4*gen-1)-P_extv).LE.cap_param%Pub_a_v)THEN
!           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-P_extv)*cap_param%alpha_v)
            radupdate=rad_v(gen)+cap_param%alpha_v*(Pressure(4*gen-1)-P_extv)*&
            (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-1)-P_extv)/(2.d0-2.d0*cap_param%num_symm_gen)
         ELSE
!            rad_update=rad_v(gen)*(1+cap_param%Pub_a_v*cap_param%alpha_v)
            radupdate=rad_v(gen)+cap_param%alpha_v*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
            /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
            /(2.d0-2.d0*cap_param%num_symm_gen)
         ENDIF
         volume_vessels=volume_vessels+&
          (2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0)

!...   Calculate Poiseuille resistance in first half of venule
           R_ven1=(8.0_dp*mu_app(gen)*L_v(gen)/2.0_dp)/(pi*radupdate**4.0_dp)

!...   CAPILLARY ELEMENT (arteriole + venule + capillary)
!...    pressure into the capillaries
           Pin_sheet=Pressure(4*gen-3)
           Pout_sheet=Pressure(4*gen-1)
!...     calulate resistance and transit time through a single capillary
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,area,area_new,alpha_c,area_scale,&
        length_scale,recruited)
            Q_sheet(gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(gen)*2**gen
            R_sheet(gen)=SHEET_RES

!...    Update soln matrix (delta p=QR across a cap)
           SparseVal(NonZeros-cap_param%num_symm_gen-1-3*(cap_param%num_symm_gen-gen))=&
            -SHEET_RES

     !...   SECOND HALF OF ARTERIOLE
!...    Update radius of arteriole based on inlet pressure NB: May need to
!....   update this later to give an average of inlet and outlet radii
      IF (Pressure(4*gen-2)-P_exta.LE.cap_param%Pub_a_v)THEN
!            radupdate=rad_a(gen)*(1+(Pressure(4*gen-2)-P_exta)*cap_param%alpha_a)
            radupdate=rad_a(gen)+cap_param%alpha_a*(Pressure(4*gen-2)-P_exta)*&
            (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-2)-P_exta)/(2.d0-2.d0*cap_param%num_symm_gen)
         ELSE
!            radupdate=rad_a(gen)*(1+cap_param%Pub_a_v*cap_param%alpha_a)
            radupdate=rad_a(gen)+cap_param%alpha_a*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
            /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
            /(2.d0-2.d0*cap_param%num_symm_gen)
         ENDIF
!         IF(nj_hypoxia.NE.0) THEN
!           IF(rad_a(gen).LE.0.25d0)THEN
!              radupdate=radupdate*k_factor
!           ENDIF
!         ENDIF
         volume_vessels=volume_vessels+&
          (2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0)

       !...  Calculate Poiseuille resistance in second half of arteriole - (only
       !...   half total generation length)
      R_art2=(8*mu_app(gen)*L_a(gen)/2.d0)/(pi*radupdate**4.d0);

!   SECOND HALF OF VENULE
!...    Update radius - linear with pressure or constant at high pressure
      IF (Pressure(4*gen)-P_extv.LE.cap_param%Pub_a_v)THEN
!           radupdate=rad_v(gen)*(1+(Pressure(4*gen)-P_extv)*cap_param%alpha_v)
            radupdate=rad_v(gen)+cap_param%alpha_v*(Pressure(4*gen)-P_extv)*&
            (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen)-P_extv)/(2.d0-2.d0*cap_param%num_symm_gen)
         ELSE
!            rad_update=rad_v(gen)*(1+cap_param%Pub_a_v*cap_param%alpha_v)
            radupdate=rad_v(gen)+cap_param%alpha_v*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
            /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
            /(2.d0-2.d0*cap_param%num_symm_gen)
         ENDIF
         volume_vessels=volume_vessels+&
          (2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0)
       !...  Poiseuille resistance in second half of venule
      R_ven2=(8.d0*mu_app(gen)*L_v(gen)/2.d0)/(pi*radupdate**4.d0)

      !...  First generation
       IF (gen.EQ.1) THEN
          SparseVal(2)=-R_art1/2.d0
          SparseVal(5)=R_ven1
          SparseVal(6)=-R_ven1/2.d0
          SparseVal(8)=-R_art2/2.d0
          SparseVal(11)=R_ven2
          SparseVal(12)=-R_ven2/2.d0
          count=12
       ELSE
        DO j=2,gen
          SparseVal(count+3+(j-2))=R_art1/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+2
      SparseVal(count)=-R_art1/2.d0**gen
        DO j=2,gen+1
          SparseVal(count+3+(j-2))=R_ven1/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+3
      SparseVal(count)=-R_ven1/2.d0**gen
         DO j=2,gen
          SparseVal(count+3+(j-2))=R_art2/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+2
      SparseVal(count)=-R_art2/2.d0**gen
        DO j=2,gen+1
          SparseVal(count+3+(j-2))=R_ven2/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+3
      SparseVal(count)=-R_ven2/2.d0**gen
      ENDIF

      ENDDO
      !...  ------------FINAL GENERATION----------------------------
      !... --------The capillaries covering the alveolar sacs---------
      !...These are just capillary beds without an associated arteriole/venule -
        Pin_sheet=Pressure(4*cap_param%num_symm_gen-6); !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*cap_param%num_symm_gen-4); !%pressure out of final capillary sheets
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,area,area_new,alpha_c,area_scale,&
        length_scale,recruited)
        SparseVal(NonZeros-cap_param%num_symm_gen-1)=-SHEET_RES
            Q_sheet(cap_param%num_symm_gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(cap_param%num_symm_gen)*2**gen
            R_sheet(cap_param%num_symm_gen)=SHEET_RES

    call enter_exit(sub_name,2)

end subroutine populate_matrix_ladder
!
!###################################################################################
!
  subroutine cap_flow_sheet(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT,&
     zone,Pin_sheet,Pout_sheet,area,area_new,alpha_c,area_scale,&
     length_scale,recruited)
!*cap_flow_sheet:* Calculates resistance across a capillary sheet.

    type(capillary_bf_parameters) :: cap_param

    integer :: ne,zone
    real(dp) :: SHEET_RES,Q_c,Hart,Hven,Pin_sheet,Pout_sheet,area,alpha_c
    real(dp) :: area_scale,length_scale

    !     Local Variables
    real(dp) :: area_new,C,Hmax_art,Hmax_ven,Hub,L_new,P_a,P_ave,&
         P_v,RBC_tt,recruited,waterfall_scale
    character(len=60) :: sub_name

    sub_name = 'cap_flow_sheet'
    call enter_exit(sub_name,1)
    !...  area of sheet is an input. Now area and pathlength need to be scaled.
    area_new=area*area_scale
    L_new=cap_param%L_c*length_scale
    !... Maximum sheet heights
    Hmax_art=cap_param%H0+alpha_c*cap_param%Pub_c
    Hmax_ven=Hmax_art
    !...  Blood pressure in and out of capillary sheet
    P_a=Pin_sheet
    P_v=Pout_sheet
    !...   CAPILLARY RECRUITMENT
    !...   CALCULATE AVERAGE CAPILLARY BLOOD PRESSURE
    P_ave=(P_a+P_v)/2.d0
    !      calculate recruitment proportion
    recruited=1.d0-cap_param%f_rec*EXP(-(P_ave**2.d0)/(cap_param%sigma_rec**2.d0))
    !      multiply flows and transit times
    area_new=area_new*recruited


    IF((P_v-cap_param%Palv).LE.cap_param%Plb_c)THEN
       ! ### Zone 2 (waterfall) scaling
       waterfall_scale=1-cap_param%F_sheet+cap_param%F_sheet* &
          EXP(-(P_v-cap_param%Palv)**2.d0/(2.d0*cap_param%sigma_cap**2.d0))
       area_new=area_new*waterfall_scale
    ENDIF
    !  -----CALCULATING SHEET FLOW----
    !...  Sheet flow model constant, C
    C=area_new/(4.d0*cap_param%mu_c*cap_param%K_cap*cap_param%f_cap*L_new**2*alpha_c)

    !###  Determine what zone we are in and calculate flow and TT
    IF((P_a-cap_param%Palv).LT.cap_param%Plb_c)THEN
       !...   ZONE 1:
       !...   Arteriole and venous pressure both less than alveolar pressure - the
       !...   sheet is collapsed
       zone=1
       Hart=0.d0
       Hven=0.d0
       Q_c=0.d0
       RBC_tt=0.d0 !TEMPORARY
    ELSEIF((P_a-cap_param%Palv).LE.cap_param%Pub_c.AND.(P_v-cap_param%Palv).LE.cap_param%Plb_c)THEN
       !...    ZONE 2:
       zone=2
       Hart=cap_param%H0+alpha_c*(P_a-cap_param%Palv)
       Hven=cap_param%H0
       Q_c=C*(Hart**4.d0-cap_param%H0**4.d0)
       RBC_tt=(3.d0*cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2.d0*alpha_c) &
            /(Hart**3.d0-cap_param%H0**3.d0)
    ELSEIF((P_a-cap_param%Palv).GT.cap_param%Pub_c.AND.(P_v-cap_param%Palv).LE.cap_param%Plb_c)THEN
       !...       ZONE 2:
       zone=2
       Hart=Hmax_art
       Hven=cap_param%H0
       Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-cap_param%Palv-cap_param%Pub_c) &
            +(Hmax_art**4.d0-cap_param%H0**4.d0)/(4.d0*alpha_c))
       RBC_tt=(3.d0*cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2.d0*alpha_c)/ &
            ((3.d0*alpha_c*Hmax_art**2.d0* &
            (P_a-cap_param%Palv-cap_param%Pub_c)+(Hmax_art**3.d0-cap_param%H0**3.d0)))
    ELSEIF((P_a-cap_param%Palv).LE.cap_param%Pub_c.AND.(P_v-cap_param%Palv).LE.cap_param%Pub_c)THEN
       !...       ZONE3 or 4: The boundary between zone 3 and 4 is not clearcut.
       zone=3 !tmp!!should = 3
       Hart=cap_param%H0+alpha_c*(P_a-cap_param%Palv)
       Hven=cap_param%H0+alpha_c*(P_v-cap_param%Palv)
       Q_c=C*(Hart**4.d0-Hven**4.d0)
       RBC_tt=(3.d0*cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2.d0*alpha_c) &
            /(Hart**3.d0-Hven**3.d0)
    ELSEIF((P_a-cap_param%Palv).GT.cap_param%Pub_c.AND.(P_v-cap_param%Palv).LE.cap_param%Pub_c)THEN
       zone=3
       Hart=Hmax_art
       Hven=cap_param%H0+alpha_c*(P_v-cap_param%Palv)
       Hub=cap_param%H0+alpha_c*cap_param%Pub_c
       Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-cap_param%Palv-cap_param%Pub_c)+ &
            (Hub**4.d0-Hven**4.d0)/(4*alpha_c))
       RBC_tt=(3.d0*cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2*alpha_c) &
            /(3.d0*alpha_c*Hmax_art**2.d0*(P_a-cap_param%Palv-cap_param%Pub_c)+ &
            (Hmax_art**3.d0-Hven**3.d0))
    ELSEIF((P_a-cap_param%Palv).GT.cap_param%Pub_c.AND.(P_v-cap_param%Palv).GT.cap_param%Pub_c)THEN
       zone=4 !!!tmp should = 3
       Hart=Hmax_art
       Hven=Hmax_ven
       Q_c=4.d0*C*alpha_c*Hart**3.d0*(P_a-P_v)
       RBC_tt=(cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
    ELSE
       !.... SOMETHING HAS GONE TERRIBLY WRONG AS THIS SHOULD BE ALL THE OPTIONS!!
       Q_c=1.d-8
       RBC_tt=0.d0
       write(*,*) 'Error, incorrect option in sheet flow model, something is wrong!', &
         P_a,P_v,cap_param%Palv,cap_param%Pub_c
    ENDIF

    !... RBC transit time (1.4 times faster than blood)
    RBC_tt=RBC_tt/1.4d0
    IF(alpha_c.LT.TOLERANCE)THEN
       C=area_new/(4.d0*cap_param%mu_c*cap_param%K_cap*cap_param%f_cap*L_new**2)
       Hart=Hmax_art
       Hven=Hmax_ven
       Q_c=4.d0*C*Hart**3.d0*(P_a-P_v)
       RBC_tt=(cap_param%mu_c*cap_param%f_cap*cap_param%K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
    ENDIF

    IF(Q_c.le.0.d0)THEN ! to acount for small negative flows when pressure in and out of a capillary are the same.
       zone=5
       Q_c=1.d-15
    ENDIF

    IF(zone.eq.1.or.zone.eq.5)THEN
       !...   ZONE 1:
       SHEET_RES=1d12
    ELSE
       !...  Resistance through a single capillary sheet, Pa.s/m^3
       SHEET_RES=(Pin_sheet-Pout_sheet)/Q_c
    ENDIF

    call enter_exit(sub_name,2)

  END subroutine cap_flow_sheet
!
!#####################
!
subroutine cap_specific_parameters(ne,Ppl,alpha_c,area_scale,length_scale,l_a,rad_a,l_v,rad_v,ngen,&
    mu_app,R_in,R_out,L_in,L_out)

    type(capillary_bf_parameters) :: cap_param

    real(dp), intent(in) :: Ppl
    integer, intent(in) :: ngen,ne
    real(dp), intent(inout) :: alpha_c,area_scale,length_scale
    real(dp), intent(inout) :: l_a(ngen),rad_a(ngen),l_v(ngen),rad_v(ngen),mu_app(ngen)
    real(dp), intent(inout) :: R_in,R_out,L_in,L_out
    character(len=60) :: sub_name

    real(dp) :: Ptp,stretch
    integer :: i

    sub_name = 'cap_specific_parameters'
    call enter_exit(sub_name,1)

!...  --CAPILLARY SHEET PROPERTIES--
!...  RELATIONSHIP BETWEEN ALVEOLAR VOLUME AND TRANSPULMONARY PRESURE
!...  Transpulmonary pressure should be defined in Pa so convert to cmH20
!...  for use in this relationship
      Ptp=(cap_param%Palv-(-Ppl))/98.06d0 !Pa -> cmH2O
!...  This assumes that zero pressure reference volume is 20% of TLC
      IF(Ptp.LT.0.d0) THEN
         length_scale=(1.0_dp-0.8_dp*exp(-0.1_dp*0.0_dp))**(1.0_dp/3.0_dp)
         area_scale=(1.0_dp-0.8_dp*exp(-0.1_dp*0.0_dp))**(2.0_dp/3.0_dp)
      ELSE
         length_scale=(1.0_dp-0.8_dp*exp(-0.1_dp*Ptp))**(1.0_dp/3.0_dp)
         area_scale=(1.0_dp-0.8_dp*exp(-0.1_dp*Ptp))**(2.0_dp/3.0_dp)
      ENDIF

!...   KSB 18/08/09: Including effect of lung inflation on sheet compliance
!...   Below value is scaled from dog lung measurements. This is valid in the
!...   human only. When we consider other animals we may need to add an
!...   option which defines alpha_c in each species.
      alpha_c=1.260_dp*((-2.04762e-09_dp*Ptp+1.3019e-07_dp)/98.06_dp) !Units m/cmH2O -> m/Pa


!    --ARTERIOLE AND VENULE PROPERTIES--
!###  APPARENT BLOOD VISCOSITY:
!...  Stepping down linearly with each generation from the larger vessels to (4e-3Pa.s) to his estimate at 30% hematocrit
!...  (1.92e-3Pa.s). (Biomechanics: Circulation)
      DO i=1,cap_param%num_symm_gen
         mu_app(i)=(4.0e-3_dp-(i-1)*(4.0e-3_dp-1.92e-3_dp)/(cap_param%num_symm_gen-1));
      ENDDO

      stretch=1.0_dp
      R_in=R_in/1000.0_dp    !radius of input artery mm->m
      R_out=R_out/1000.0_dp  ! radius of outlet vein mm->m
      L_in=L_in/1000.0_dp
      L_out=L_out/1000.0_dp
!...  Note that stretch is included here so we don't have to add it later
      do i=1,cap_param%num_symm_gen
         L_a(i)=L_in-i*(L_in-cap_param%L_art_terminal*stretch)/cap_param%num_symm_gen
         L_v(i)=L_out-i*(L_in-cap_param%L_vein_terminal*stretch)/cap_param%num_symm_gen
      enddo
      do i=1,cap_param%num_symm_gen
        rad_a(i)=(cap_param%R_art_terminal*sqrt(1.0_dp/stretch)-R_in)*i/cap_param%num_symm_gen+R_in
        rad_v(i)=(cap_param%R_vein_terminal*sqrt(1.0_dp/stretch)-R_out)*i/cap_param%num_symm_gen+R_out
      ENDDO
      cap_param%alpha_a=R_in/(6670.0_dp)
      cap_param%alpha_v=R_out/(6670.0_dp)
      !store element length
      elem_field(ne_length,ne)=L_a(1)*1000.0_dp

    call enter_exit(sub_name,2)
end subroutine cap_specific_parameters
!
!################################################
!
subroutine cap_flow_admit(ne,admit,eff_admit_downstream,Lin,Lout,P1,P2,&
  Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,no_freq,harmonic_scale,elast_param,&
  cap_model)

  integer,intent(in) :: ne
  integer,intent(in) :: no_freq
  complex(dp), intent(inout) :: admit(no_freq)
  complex(dp), intent(in) :: eff_admit_downstream(no_freq)
  real(dp), intent(inout) ::Lin,Lout,P1,P2,Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap
  real(dp), intent(in) :: harmonic_scale
  integer, intent(in) :: cap_model

  type(elasticity_param) :: elast_param

  type(capillary_bf_parameters) :: cap_param
  type(fluid_properties) :: fp
  integer :: ngen
  real(dp) :: alpha_c,area_scale,length_scale
  real(dp) :: radupdate,P_exta,P_extv,R_art1,R_ven1,R_art2,R_ven2,Q01_mthrees,Pin,Pout
  integer :: gen
  real(dp) :: SHEET_RES,Q_c,Hart,Hven,RBC_TT,area,area_new,recruited,omega,nd_omega
  integer :: zone,nf
  complex(dp) :: bessel0,bessel1,f10
  real(dp) :: wolmer,wavespeed
  integer :: i,iter,j,num_sheet
  integer, allocatable :: SparseCol(:)
  integer, allocatable ::SparseRow(:)
  real(dp) :: ErrorEstimate,Pin_SHEET,Pout_SHEET
  real(dp), allocatable :: Pressure(:)
  real(dp) ::  Qtot,Qgen
  real(dp),allocatable :: Q_sheet(:)
  real(dp),allocatable :: RHS(:)
  real(dp) :: Rtot
  real(dp),allocatable :: Solution(:)
  real(dp),allocatable :: SolutionLast(:)
  real(dp),allocatable :: SparseVal(:)
  real(dp) :: TOTAL_CAP_VOL,TOTAL_SHEET_H,TOTAL_SHEET_SA,&
         TT_TOTAL,R_upstream,R_downstream
  real(dp),allocatable :: l_a(:),rad_a(:),l_v(:),rad_v(:),mu_app(:)
  integer :: MatrixSize,NonZeros,submatrixsize
  real(dp) :: sheet_number
  integer SOLVER_FLAG

  complex(dp), allocatable :: cap_admit(:,:)
  complex(dp), allocatable :: tube_admit(:,:)
  complex(dp), allocatable :: cap_eff_admit(:,:)
  complex(dp), allocatable :: tube_eff_admit(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp), allocatable :: prop_const_cap(:,:)
  complex(dp), allocatable :: sheet_admit_matrix(:,:,:)

  complex(dp) :: daughter_admit,sister_admit,reflect_coeff,sister_current


  character(len=60) :: sub_name
  integer :: AllocateStatus

  sub_name = 'cap_flow_admit'
  call enter_exit(sub_name,1)

  !     Number of non-zero entries in solution matrix.
  NonZeros=3
  do i=2,cap_param%num_symm_gen
      NonZeros=NonZeros+4*i+10
  enddo
  !     The size of the solution matrix (number of unknown pressures and flows)
  MatrixSize=5*cap_param%num_symm_gen-3
  !!     The number of unknown pressures
  submatrixsize=4*cap_param%num_symm_gen-4
  ngen=cap_param%num_symm_gen

  !...  ---INITIALISATION
  !...  The input Q01 gives us an estimate for flow into the acinus from the large
  !...  vessel model.
  !...  This is in mm^3/s and needs to be converted to m^3/s to use in calculating
  !...  arteriole and venule resistance
  Q01_mthrees=Q01/1.d9 !mm3/s->m3/s
  Pin = P1
  Pout = P2
  !     Sheet area (unscaled):
  !...  We define a sheet area for input into the capillary model.
  !...  This area is at full inflation and will be scaled within CAP_FLOW_SHEET
  !...  Area of an individual sheet
  sheet_number=0
  DO i=1,cap_param%num_symm_gen
      sheet_number=sheet_number+2.d0**i
  ENDDO
  area=cap_param%total_cap_area/(sheet_number*num_units) !m^2
  ngen=cap_param%num_symm_gen
  !1. need to resolve the micro-unit flow to get correct pressure and flow in each unit
  allocate (Pressure(submatrixsize), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array Press ***"
  allocate (SparseCol(NonZeros), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array SparseCol***"
  allocate (SparseVal(NonZeros), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array  SParseVal***"
  allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array SparseRow ***"
  allocate (Solution(MatrixSize), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array Solution ***"
  allocate (SolutionLast(MatrixSize), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array SolutionLast***"
  allocate (RHS(MatrixSize), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array RHS***"
  allocate (l_a(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array l_a***"
  allocate (rad_a(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array rad_a***"
  allocate (l_v(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array l_v ***"
  allocate (rad_v(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array rad_v ***"
  allocate (Q_Sheet(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array Q_sheet***"
  allocate (mu_app(ngen), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array mu_app***"

  Pressure=0.0_dp
  SparseCol=0
  SparseRow=0
  SparseVal=0.0_dp
  Solution=0.0_dp
  SolutionLast=0.0_dp

  !!...  Initial guess for pressure distribution lets say all arterial pressures are the same
  !!...  and all the venous pressures are the same solution appears independent of this.
  DO i=1,cap_param%num_symm_gen-1
      Pressure(4*i-3)=1000.0_dp ! Pa
      Pressure(4*i-2)=1000.0_dp
      Pressure(4*i-1)=100.0_dp
      Pressure(4*i)=100.0_dp
  ENDDO


  !###  INITIAL SOLUTION GUESS
  DO i=1,submatrixsize
      solution(i)=Pressure(i)
  ENDDO
  DO i=submatrixsize+1,matrixSize-1
      solution(i)=Q01_mthrees/2**cap_param%num_symm_gen
  ENDDO
  Solution(Matrixsize)=Q01_mthrees
  !###  INITIALISE SOLUTIONLAST
  DO j=1,MatrixSize
      SolutionLast(j)=Solution(j)
  ENDDO

  !### INPUT TO THE LADDER MODEL THAT IS INDEPENDENT OF ITERATION
  CALL LADDERSOL_MATRIX(NonZeros,MatrixSize,submatrixsize,&
      SparseCol,SparseRow,SparseVal,RHS,Pin,Pout)

  CALL cap_specific_parameters(ne,Ppl,alpha_c,area_scale,length_scale,l_a,rad_a,l_v,rad_v,ngen,&
      mu_app,Rin,Rout,Lin,Lout)
  !### ITERATIVE LOOP
  iter=0
  ErrorEstimate=1.d10
  DO WHILE(ErrorEstimate.GT.1.0d-9.AND.iter.LT.100)
      iter=iter+1
      !...  CALCULATE RESISTANCE GIVEN CURRENT PRESSURE AND FLOW - THEN UPDATE
      !...  SparseVal- THese are the only elements of the solution matrix that need
      !.... iteratively updating

      CALL POPULATE_MATRIX_LADDER(ne,NonZeros,submatrixsize,ngen,area,alpha_c,&
          area_scale,length_scale,mu_app,Ppl,Pin,Pout,Pressure,&
          Q01_mthrees,Q_sheet,SparseVal,l_a,rad_a,l_v,rad_v)


      call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
          Solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)


      DO j=1,submatrixsize
          Pressure(j)=Solution(j)
      ENDDO
      Q01_mthrees=Solution(MatrixSize)
      !     Estimating Error in solution
      ErrorEstimate=0.d0
      DO i=1,MatrixSize
          ErrorEstimate=ErrorEstimate+&
              DABS((Solution(i)-SolutionLast(i))**2.d0&
              /Solution(i)**2.d0)
          SolutionLast(i)=Solution(i)
      ENDDO
      ErrorEstimate=ErrorEstimate/MatrixSize

  ENDDO

  Qtot=0
  Do i=1,cap_param%num_symm_gen
     Qtot=Qtot+Q_sheet(i)*2.d0**i
  ENDDO

  deallocate (SparseCol, STAT = AllocateStatus)
  deallocate (SparseRow, STAT = AllocateStatus)
  deallocate (SparseRow, STAT = AllocateStatus)
  deallocate (Solution, STAT = AllocateStatus)
  deallocate (SolutionLast, STAT = AllocateStatus)
  deallocate (RHS, STAT = AllocateStatus)

  allocate (cap_admit(ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for cap_admit ***"
  allocate(tube_admit(4*ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for tube_admit ***"
  allocate (cap_eff_admit(ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for cap_eff_admit ***"
  allocate(tube_eff_admit(4*ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for tube_eff_admit ***"
  allocate(prop_const(4*ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const ***"
  allocate(prop_const_cap(ngen,no_freq), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const_cap ***"
  if(cap_model == 2)then !Allocating memory for sheet admittance matrix
    allocate(sheet_admit_matrix(ngen,no_freq,4), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const_cap ***"
  endif


  cap_admit=cmplx(0.0_dp,0.0_dp,8) !initialise
  tube_admit=cmplx(0.0_dp,0.0_dp,8) !initialise
  cap_eff_admit=cmplx(0.0_dp,0.0_dp,8) !initialise
  tube_eff_admit=cmplx(0.0_dp,0.0_dp,8) !initialise
  prop_const=cmplx(0.0_dp,0.0_dp,8) !initialise
  prop_const_cap=cmplx(0.0_dp,0.0_dp,8) !initialise


  radupdate=0.d0
  do gen=1,cap_param%num_symm_gen-1
    Q01_mthrees=Qtot/2.0_dp
!!...    FIRST HALF OF ARTERIOLE
    if(rad_a(gen).LT.100.d-6) then
      P_exta=cap_param%Palv ! From Yen Alveolar pressure dominates vessels <200um diam
    else
      P_exta=-Ppl
    endif
    if ((Pin-P_exta).LE.cap_param%Pub_a_v)THEN
      radupdate=rad_a(gen)+cap_param%alpha_a*(Pin-P_exta)*&
         (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
         (Pin-P_exta)/(2.d0-2.d0*cap_param%num_symm_gen)
    else
      radupdate=rad_a(gen)+cap_param%alpha_a*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
        /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
        /(2.d0-2.d0*cap_param%num_symm_gen)
    endif
!!...  Calculate Poiseuille resistance in first half of arteriole - (only
!!...  half total generation length)
     R_art1=(8.d0*mu_app(gen)*L_a(gen)/2.d0)/(PI*radupdate**4.d0)
     Pin=Pin-R_art1*Q01_mthrees
     do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
      !!!ARC TO FIX alpha_a is in m/Pa, need in 1/Pa (just read in from main model?)
       omega=nf*2*PI*harmonic_scale
       wolmer=(radupdate*1000.0_dp)*sqrt(omega*fp%blood_density/mu_app(gen))
       call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
       f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)
       wavespeed=sqrt(1.0_dp/(2*fp%blood_density*(elast_param%elasticity_parameters(1))))*sqrt(1-f10)!alpha in the sense of this model is 1/Pa so has to be dovided by radius
       tube_admit(gen,nf)=PI*(radupdate*1000.0_dp)**2/(fp%blood_density*wavespeed/sqrt(1-f10))*sqrt(1-f10)
       prop_const(gen,nf)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed)
     enddo
!!...    FIRST HALF OF VENULE
!!...    Update radius of venule based on outlet pressure
     if(rad_v(gen).LT.100.d-6) THEN
       P_extv=cap_param%Palv ! From Yen Alveolar pressure dominates vessels <200um diam
     else
       P_extv=-Ppl
     endif
     if((Pout-P_extv).LE.cap_param%Pub_a_v)THEN
       radupdate=rad_v(gen)+cap_param%alpha_v*(Pout-P_extv)*&
         (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
         (Pout-P_extv)/(2.d0-2.d0*cap_param%num_symm_gen)
     else
       radupdate=rad_v(gen)+cap_param%alpha_v*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
         /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
         /(2.d0-2.d0*cap_param%num_symm_gen)
     endif
!!...   Calculate Poiseuille resistance in first half of venule
     R_ven1=(8.0_dp*mu_app(gen)*L_v(gen)/2.0_dp)/(pi*radupdate**4.0_dp)
     Pout=Pout+R_ven1*Q01_mthrees
     do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
       omega=nf*2*PI*harmonic_scale
       wolmer=(radupdate*1000.0_dp)*sqrt(omega*fp%blood_density/mu_app(gen))
       call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
       f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)
       wavespeed=sqrt(1.0_dp/(2*fp%blood_density*elast_param%elasticity_parameters(1)))*sqrt(1-f10) !mm/s
       tube_admit(gen+2*ngen,nf)=PI*(radupdate*1000.0_dp)**2/(fp%blood_density*wavespeed/sqrt(1-f10))*sqrt(1-f10)!mm3/Pa.s
       prop_const(gen+2*ngen,nf)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed)!1/mm
     enddo

!!...   CAPILLARY ELEMENT (arteriole + venule + capillary)
     Q01_mthrees=Q01_mthrees-Q_sheet(gen)
     Hart=cap_param%H0+alpha_c*(Pressure(4*gen-3)-cap_param%Palv)
     Hven=cap_param%H0+alpha_c*(Pressure(4*gen-1)-cap_param%Palv)
    do nf = 1,no_freq
      omega=nf*2*PI*harmonic_scale
      if(cap_model == 1)then!
        call calc_cap_admit_consth(Hart,Hven,omega,alpha_c,cap_admit(gen,nf),&
          prop_const_cap(gen,nf))
      elseif(cap_model == 2)then
        call calc_cap_admit_varh(Hart,Hven,omega,alpha_c,&
          sheet_admit_matrix(gen,nf,1),sheet_admit_matrix(gen,nf,2),&
          sheet_admit_matrix(gen,nf,3),sheet_admit_matrix(gen,nf,4), &
          prop_const_cap(gen,nf)) ! Non-constant capillary sheet
      endif
    enddo
     !...   SECOND HALF OF ARTERIOLE
!...    Update radius of arteriole based on inlet pressure
     if(Pin-P_exta.LE.cap_param%Pub_a_v)THEN
       radupdate=rad_a(gen)+cap_param%alpha_a*(Pin-P_exta)*&
         (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
         (Pin-P_exta)/(2.d0-2.d0*cap_param%num_symm_gen)
     else
       radupdate=rad_a(gen)+cap_param%alpha_a*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
         /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
         /(2.d0-2.d0*cap_param%num_symm_gen)
     endif

       !...  Calculate Poiseuille resistance in second half of arteriole - (only
       !...   half total generation length)
      R_art2=(8*mu_app(gen)*L_a(gen)/2.d0)/(pi*radupdate**4.d0)
      Pin=Pin-R_art2*Q01_mthrees

     do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
       omega=nf*2*PI*harmonic_scale
       wolmer=(radupdate*1000.0_dp)*sqrt(omega*fp%blood_density/mu_app(gen))
       call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
       f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)
       wavespeed=sqrt(1.0_dp/(2*fp%blood_density*elast_param%elasticity_parameters(1)))*sqrt(1-f10)
       tube_admit(gen+ngen,nf)=PI*(radupdate*1000.0_dp)**2/(fp%blood_density*wavespeed/sqrt(1-f10))*sqrt(1-f10)
       prop_const(gen+ngen,nf)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed)
     enddo

!   SECOND HALF OF VENULE
!...    Update radius - linear with pressure or constant at high pressure
     if (Pout-P_extv.LE.cap_param%Pub_a_v)THEN
       radupdate=rad_v(gen)+cap_param%alpha_v*(Pout-P_extv)*&
         (gen-cap_param%num_symm_gen)/(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*&
          (Pout-P_extv)/(2.d0-2.d0*cap_param%num_symm_gen)
     else
      radupdate=rad_v(gen)+cap_param%alpha_v*cap_param%Pub_a_v*(gen-cap_param%num_symm_gen)&
        /(1.d0-cap_param%num_symm_gen)+alpha_c*(1-gen)*cap_param%Pub_a_v&
        /(2.d0-2.d0*cap_param%num_symm_gen)
     endif
     R_ven2=(8.d0*mu_app(gen)*L_v(gen)/2.d0)/(pi*radupdate**4.d0)
     Pout=Pout-R_ven2*Q01_mthrees
     do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
       omega=nf*2*PI*harmonic_scale
       wolmer=(radupdate*1000.0_dp)*sqrt(omega*fp%blood_density/mu_app(gen))
       call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
       f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)
       wavespeed=sqrt(1.0_dp/(2*fp%blood_density*elast_param%elasticity_parameters(1)))*sqrt(1-f10) !mm/s
       tube_admit(gen+3*ngen,nf)=PI*(radupdate*1000.0_dp)**2/(fp%blood_density*wavespeed/sqrt(1-f10))*sqrt(1-f10)
       prop_const(gen+3*ngen,nf)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed) !1/mm
     enddo
  enddo
  !!...   CAPILLARY ELEMENT (arteriole + venule + capillary)  at terminal
    Pin_sheet=Pressure(4*cap_param%num_symm_gen-6) !pressure into final capillary sheets
    Pout_sheet=Pressure(4*cap_param%num_symm_gen-4) !pressure out of final capillary sheets
    Hart=cap_param%H0+alpha_c*(Pin_sheet-cap_param%Palv) !Sheet height at arterial side
    Hven=cap_param%H0+alpha_c*(Pout_sheet-cap_param%Palv) !sheet height at venous side
    gen=cap_param%num_symm_gen
    do nf = 1,no_freq
      omega=nf*2*PI*harmonic_scale
      if(cap_model == 1)then!
        call calc_cap_admit_consth(Hart,Hven,omega,alpha_c,cap_admit(gen,nf),&
          prop_const_cap(gen,nf))
      elseif(cap_model == 2)then
        call calc_cap_admit_varh(Hart,Hven,omega,alpha_c,&
          sheet_admit_matrix(gen,nf,1),sheet_admit_matrix(gen,nf,2),&
          sheet_admit_matrix(gen,nf,3),sheet_admit_matrix(gen,nf,4), &
          prop_const_cap(gen,nf)) ! Non-constant capillary sheet
      endif
    enddo

  !FIRST GENERATION - need to relate to vein outlet admittance
  !first vein in generation
  do nf=1,no_freq
  !sister and vessel are identical H=1.0
   reflect_coeff=(2*tube_admit(1+2*ngen,nf)-eff_admit_downstream(nf))/(eff_admit_downstream(nf)+2*tube_admit(1+2*ngen,nf))
   tube_eff_admit(1+2*ngen,nf)=tube_admit(1+2*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(1+2*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(1+2*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))
  enddo
  !Now second vein and then capillary
  do nf=1,no_freq
   !Vein and capillary are sisters
    if(cap_model == 1)then!
      sister_current=exp(-1.0_dp*prop_const_cap(1,nf)*cap_param%L_c*1000.0_dp)/&
        exp(-1.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp)
      reflect_coeff=(tube_admit(1+3*ngen,nf)+(2*sister_current-1)*cap_admit(1,nf)-tube_admit(1+2*ngen,nf))&
        /(tube_admit(1+2*ngen,nf)+tube_admit(1+3*ngen,nf)+cap_admit(1,nf))
      tube_eff_admit(1+3*ngen,nf)=tube_admit(1+3*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))
      !sister is the tube current is the capillary
      sister_current=exp(-1.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp)/&
        exp(-1.0_dp*prop_const_cap(1,nf)*cap_param%L_c*1000.0_dp)
      reflect_coeff=((2*sister_current-1)*tube_admit(1+3*ngen,nf)+cap_admit(1,nf)-tube_admit(1+2*ngen,nf))&
        /(tube_admit(1+2*ngen,nf)+tube_admit(1+3*ngen,nf)+cap_admit(1,nf))
      cap_eff_admit(1,nf)=cap_admit(1,nf)*(1&
        -reflect_coeff*exp(-2.0_dp*prop_const_cap(1,nf)*cap_param%L_c*1000.0_dp))/&
        (1+reflect_coeff*exp(-2.0_dp*prop_const_cap(1,nf)*cap_param%L_c*1000.0_dp))
    elseif(cap_model == 2)then
      ! Take the capillary first
      !Ycap = Y11-Y12Y21/(Ydaughter-Ysister+Y22)
      cap_eff_admit(1,nf) = sheet_admit_matrix(1,nf,1)&
          -(sheet_admit_matrix(1,nf,2)*sheet_admit_matrix(1,nf,3))/&
          (tube_admit(1+2*ngen,nf)- tube_admit(1+3*ngen,nf)+sheet_admit_matrix(1,nf,4))
      ! Now the vein
      !Calculate the capillary admittance at the venous side
      sister_admit = sheet_admit_matrix(1,nf,3)*sheet_admit_matrix(1,nf,2)/ &
        (cap_eff_admit(1,nf) - sheet_admit_matrix(1,nf,3)) + sheet_admit_matrix(1,nf,4)
      !! Calculate the reflection coefficient
      sister_current = exp(-1.0_dp*prop_const_cap(1,nf)*cap_param%L_c*1000.0_dp)/&
        exp(-1.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp)
      reflect_coeff=(tube_admit(1+3*ngen,nf)+(2*sister_current-1)*sister_admit-tube_admit(1+2*ngen,nf))&
        /(tube_admit(1+2*ngen,nf)+tube_admit(1+3*ngen,nf)+sister_admit)
      tube_eff_admit(1+3*ngen,nf)=tube_admit(1+3*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(1+3*ngen,nf)*L_v(1)*1000.0_dp/2.0_dp))
    endif
  enddo

  do gen=2,cap_param%num_symm_gen-1
    do nf =1,no_freq
      !Next two vessels are identical
      reflect_coeff=(2*tube_admit(gen+2*ngen,nf)-tube_admit(gen+2*ngen-1,nf))&
        /(tube_admit(gen+2*ngen-1,nf)+2*tube_admit(gen+2*ngen,nf))
       tube_eff_admit(gen+2*ngen,nf)=tube_admit(gen+2*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(gen+2*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(gen+2*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))
      if(cap_model == 1)then!
        !Next up a vein plus a capillary
        sister_current=exp(-1.0_dp*prop_const_cap(gen,nf)*cap_param%L_c*1000.0_dp)/&
          exp(-1.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp)
        reflect_coeff=(tube_admit(gen+3*ngen,nf)+(2*sister_current-1)*cap_admit(gen,nf)-tube_admit(gen+2*ngen,nf))&
          /(tube_admit(gen+2*ngen,nf)+tube_admit(gen+3*ngen,nf)+cap_admit(gen,nf))
        tube_eff_admit(gen+3*ngen,nf)=tube_admit(gen+3*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))
        !sister is the tube current is the capillary
        sister_current=exp(-1.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp)/&
          exp(-1.0_dp*prop_const_cap(gen,nf)*cap_param%L_c*1000.0_dp)
        reflect_coeff=((2*sister_current-1)*tube_admit(gen+3*ngen,nf)+cap_admit(gen,nf)-tube_admit(1+2*ngen,nf))&
           /(tube_admit(1+2*ngen,nf)+tube_admit(gen+3*ngen,nf)+cap_admit(gen,nf))
        cap_eff_admit(gen,nf)=cap_admit(gen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const_cap(gen,nf)*cap_param%L_c*1000.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const_cap(gen,nf)*cap_param%L_c*1000.0_dp))
      elseif(cap_model == 2)then
        ! Take the capillary first
        !Ycap = Y11-Y12Y21/(Ydaughter-Ysister+Y22)
        cap_eff_admit(gen,nf) = sheet_admit_matrix(gen,nf,1)&
          -(sheet_admit_matrix(gen,nf,2)*sheet_admit_matrix(gen,nf,3))/&
          (tube_admit(gen+2*ngen,nf)- tube_admit(gen+3*ngen,nf)+sheet_admit_matrix(gen,nf,4))
        ! Now the vein
        !Calculate the capillary admittance at the venous side
        sister_admit = sheet_admit_matrix(gen,nf,3)*sheet_admit_matrix(gen,nf,2)/ &
          (cap_eff_admit(gen,nf) - sheet_admit_matrix(gen,nf,3)) + sheet_admit_matrix(gen,nf,4)
        !! Calculate the reflection coefficient
        sister_current=exp(-1.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp)/&
          exp(-1.0_dp*prop_const_cap(gen,nf)*cap_param%L_c*1000.0_dp)
        reflect_coeff=(tube_admit(gen+3*ngen,nf)+(2*sister_current-1)*sister_admit-tube_admit(gen+2*ngen,nf))&
          /(tube_admit(gen+2*ngen,nf)+tube_admit(gen+3*ngen,nf)+sister_admit)
        tube_eff_admit(gen+3*ngen,nf)=tube_admit(gen+3*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(gen+3*ngen,nf)*L_v(gen)*1000.0_dp/2.0_dp))
      endif
    enddo
  enddo

 !Final generation capillary
   gen=cap_param%num_symm_gen-1
   do nf=1,no_freq
     if(cap_model == 1)then!
       !two identical sisters
       reflect_coeff=(2*cap_admit(cap_param%num_symm_gen,nf)-tube_admit(gen+3*ngen,nf))&
         /(tube_admit(gen+3*ngen,nf)+2*cap_admit(cap_param%num_symm_gen,nf))
       cap_eff_admit(cap_param%num_symm_gen,nf)=cap_admit(cap_param%num_symm_gen,nf)*(1&
         -reflect_coeff*exp(-2.0_dp*prop_const_cap(cap_param%num_symm_gen,nf)**cap_param%L_c*1000.0_dp))/&
         (1+reflect_coeff*exp(-2.0_dp*prop_const_cap(cap_param%num_symm_gen,nf)*cap_param%L_c*1000.0_dp))
     elseif(cap_model == 2)then
        !two identical sisters
        !Ycap = Y11-2Y12Y21/(Ytube+2Y22)
        cap_eff_admit(cap_param%num_symm_gen,nf) = sheet_admit_matrix(cap_param%num_symm_gen,nf,1)&
          -(2.0_dp*sheet_admit_matrix(cap_param%num_symm_gen,nf,2)*sheet_admit_matrix(cap_param%num_symm_gen,nf,3))/&
          (tube_admit(gen+3*ngen,nf)+2.0_dp*sheet_admit_matrix(cap_param%num_symm_gen,nf,4))
     endif
  enddo

  !now calculate effective admittance up the arteriole side of the tree
   do gen=cap_param%num_symm_gen-1,1,-1
     do nf=1,no_freq
   !Second half of ateriole daughter admittance is 2* vessel below it
     if(gen.eq.(cap_param%num_symm_gen-1))then
       daughter_admit=2.0_dp*cap_eff_admit(cap_param%num_symm_gen,nf)
     else
       daughter_admit=2.0_dp*tube_eff_admit((gen+1)+2*ngen,nf)
     endif
     reflect_coeff=(tube_admit(gen+2*ngen,nf)-daughter_admit)/&
              (tube_admit(gen+2*gen,nf)+daughter_admit)
     tube_eff_admit(gen+2*ngen,nf)=tube_admit(gen+2*ngen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(gen+2*ngen,nf)*L_a(gen)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(gen+2*ngen,nf)*L_a(gen)*1000.0_dp/2.0_dp))
      !first half of arteriole, daughter admittance is second half of arteriole plus cap
      daughter_admit=tube_eff_admit(gen+2*ngen,nf)+cap_eff_admit(gen,nf)
      reflect_coeff=(tube_admit(gen,nf)-daughter_admit)/&
              (tube_admit(gen,nf)+daughter_admit)
      tube_eff_admit(gen,nf)=tube_admit(gen,nf)*(1&
              -reflect_coeff*exp(-2.0_dp*prop_const(gen,nf)*L_a(gen)*1000.0_dp/2.0_dp))/&
              (1+reflect_coeff*exp(-2.0_dp*prop_const(gen,nf)*L_a(gen)*1000.0_dp/2.0_dp))

    enddo
   enddo

   do nf=1,no_freq
     admit(nf)=tube_eff_admit(1,nf)!tube_eff_admit(1,nf) !!effective daughter admittance of first arteriole
   enddo

  deallocate(cap_admit)
  deallocate(tube_admit)
  deallocate(cap_eff_admit)
  deallocate(tube_eff_admit)
  deallocate(prop_const)
  if(cap_model == 2)then
    deallocate(sheet_admit_matrix)
  endif
  deallocate (Pressure, STAT = AllocateStatus)
  deallocate (l_a, STAT = AllocateStatus)
  deallocate (l_v, STAT = AllocateStatus)
  deallocate (rad_a, STAT = AllocateStatus)
  deallocate (rad_v, STAT = AllocateStatus)
  deallocate (Q_sheet, STAT = AllocateStatus)
  deallocate (mu_app, STAT = AllocateStatus)

  !write(*,*) 'Capillary impedance completed for element: ',ne
  call enter_exit(sub_name,2)

end subroutine cap_flow_admit

!
!################################################################
!
subroutine calc_cap_admit_consth(Hart,Hven,omega,alpha_c,Y,prop_const)

  ! Parameters:
  real(dp), intent(in) :: Hart,Hven,omega,alpha_c
  complex(dp), intent(out) :: Y, prop_const

  !Local variables
  complex(dp) :: Gamma_sheet
  type(capillary_bf_parameters) :: cap_param

  character(len=60) :: sub_name
  sub_name = 'calc_cap_admit_consth'
  call enter_exit(sub_name,1)

    Gamma_sheet=sqrt(omega*Hart**3*cap_param%L_c**2*alpha_c/&
        (cap_param%mu_c*cap_param%K_cap*cap_param%F_cap))*&
        sqrt(cmplx(0.0_dp,1.0_dp,8))*1000.0_dp**3
    Y=Gamma_sheet
    prop_const=sqrt(cmplx(0.0_dp,1.0_dp,8)*cap_param%mu_c* &
    cap_param%K_cap*cap_param%F_cap*omega*alpha_c/(Hart**3))/1000.0_dp!1/mm

  call enter_exit(sub_name,2)

end subroutine calc_cap_admit_consth
!
!################################################################
!
subroutine calc_cap_admit_varh(Hart,Hven,omega_d,alpha_c,Y11,Y12,Y21,Y22,prop_const)
! Calculating the admittance components for non-constant sheet height
! Used by cap_flow_admit to solve the time dependent capillary sheet impedance.
!    Inputs:
!           1- ha: Dimensional sheet height at arterial side of the capillary
!           2- hv: Dimensional sheet height at venous side of the capillary
!           3- omega: Dimensionless oscillatory frequency
!    Output:
!           1- As of now it will print out the constants from Fung's paper (Pulmonary microvascular impedance-1972)
!           2- Two port network capillary input admittance components


! Parameters:
real(dp), intent(in) :: Hart,Hven,omega_d,alpha_c
complex(dp), intent(out) :: Y11, Y12, Y21, Y22,prop_const


! Local Variables:
type(capillary_bf_parameters) :: cap_param
integer, parameter :: N_nodes = 2500
integer :: n,ldb
integer :: iopt, info
real(dp),  allocatable :: sparseval(:)
integer, allocatable :: sparsecol(:), sparserow(:)
real(dp) :: stiff1(2*N_nodes+2,2*N_nodes+2), stiff2(2*N_nodes+2,2*N_nodes+2)
real(dp) :: RHS1(2*N_nodes+2), RHS2(2*N_nodes+2)
integer :: nrhs
real(dp) :: deriv1(2),deriv2(2),dx,h(5),rep
integer :: i
integer :: factors(8)
integer :: NonZeros
complex(dp) :: C1,C2,C3,C4
real(dp) :: ha,hv,omega !Non-dimensional parameters
character(len=60) :: sub_name

sub_name = 'calc_cap_admit_varh'
call enter_exit(sub_name,1)

prop_const = 0.0
Y11 = 0.0
Y12 = 0.0
Y21 = 0.0
Y22 = 0.0
#ifdef HAVE_SUPERLU
!Convert dimensional parameters to non-dimensional ones
omega = (cap_param%mu_c*cap_param%K_cap*cap_param%F_cap*omega_d*alpha_c*(cap_param%L_c)**2)/(cap_param%H0**3) ! Non-Dimensional Frequency
ha = Hart/cap_param%H0
hv = Hven/cap_param%H0

rep = ha**4 - hv**4
dx = 1.0/(N_nodes - 1.0)
h(1) = ha**3
h(2) = (ha**4 - 2*rep*dx)**(0.75)
h(3) = (ha**4 - (N_nodes-1)*rep*dx)**(0.75)
h(4) = (ha**4 - (N_nodes-2)*rep*dx)**(0.75)
h(5) = (ha**4 - (N_nodes-3)*rep*dx)**(0.75)
Y11 = 0 ! admittance initialisation
Y12 = 0 ! admittance initialisation
Y21 = 0 ! admittance initialisation
Y22 = 0 ! admittance initialisation
n = N_Nodes*2 + 2
call Matrix(N_nodes, ha, hv, omega, stiff1, stiff2, RHS1, RHS2)
call Mat_to_CC(stiff1,N_nodes,sparsecol,sparserow,sparseval,NonZeros)
 nrhs = 1
  ldb = n

!
!  Factor the matrix.
!
  iopt = 1
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS1, ldb, factors, info )

!
!  Solve the factored system.
!
  iopt = 2
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS1, ldb, factors, info )

!write (*,*) 'C_1 = ',RHS1(2),'+ i(',RHS1(N_Nodes+2),')'
 C1 = RHS1(2) + RHS1(N_Nodes+2)*cmplx(0.0_dp,1.0_dp,8)
deriv1(1) = (h(2)*RHS1(3) - h(1)*RHS1(1))/(2*dx)
deriv1(2) = (h(2)*RHS1(N_nodes+4) - h(1)*RHS1(N_Nodes+2))/(2*dx)
!write (*,*) 'C_2 = ',deriv1(1),'+ i(',deriv1(2),')'
 C2 = deriv1(1) + deriv1(2)*cmplx(0.0_dp,1.0_dp,8)

!
!  Free memory.
!
  iopt = 3
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS1, ldb, factors, info )
!
!  Terminate.
!

call Mat_to_CC(stiff2,N_nodes,sparsecol,sparserow,sparseval,NonZeros)
 nrhs = 1
  ldb = n

!
!  Factor the matrix.
!
  iopt = 1
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS2, ldb, factors, info )

!
!  Solve the factored system.
!
  iopt = 2
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS2, ldb, factors, info )

!write (*,*) 'C_3 = ',RHS2(N_Nodes+1),'+ i(',RHS2(2*N_Nodes+2),')'
 C3 = RHS2(N_Nodes+1) + RHS2(2*N_Nodes+2)*cmplx(0.0_dp,1.0_dp,8)
deriv2(1) = (3*h(3)*RHS2(N_nodes+1) - 4*h(4)*RHS2(N_nodes) + h(5)*RHS2(N_nodes-1))/(2*dx)
deriv2(2) = (3*h(3)*RHS2(2*N_nodes+2) - 4*h(4)*RHS2(2*N_nodes+1) + h(5)*RHS2(2*N_nodes))/(2*dx)
!write (*,*) 'C_4 = ',deriv2(1),'+ i(',deriv2(2),')'
 C4 = deriv2(1) + deriv2(2)*cmplx(0.0_dp,1.0_dp,8)

!
!  Free memory.
!
  iopt = 3
  call c_fortran_dgssv ( iopt, n, NonZeros, nrhs, sparseval, sparserow, &
    sparsecol, RHS1, ldb, factors, info )
!
!  Terminate.
!

Y11 = -C2/C1
Y12 = 1.0/C3
Y21 = 1.0/C1
Y22 = -C4/C3
!write(*,*) '=================================================='
!write(*,*) 'C1: ', C1*C2
!write(*,*) 'IMPEDANCE: ', Y11 - (Y12*Y21)/(Y22)

!Y11, Y12, Y21, Y22 are all dimensionless, but for the rest of the model they should have units of admittance
!Here we need to convert back to dimensional form
Y11 = Y11*cap_param%H0**3/(cap_param%mu_c*cap_param%K_cap)*1000.0_dp**3 !m->mm
Y12 = Y12*cap_param%H0**3/(cap_param%mu_c*cap_param%K_cap)*1000.0_dp**3 !m->mm
Y21 = Y21*cap_param%H0**3/(cap_param%mu_c*cap_param%K_cap)*1000.0_dp**3 !m->mm
Y22 = Y22*cap_param%H0**3/(cap_param%mu_c*cap_param%K_cap)*1000.0_dp**3 !m->mm

prop_const=sqrt(cmplx(0.0_dp,1.0_dp,8)*cap_param%mu_c* &
    cap_param%K_cap*cap_param%F_cap*omega*alpha_c/ &
    (((Hart+Hven)/2.0_dp)**3.0_dp))/1000.0_dp!1/mm

#else
write(*,*) 'SuperLU is not available you cannot use capillary model 2.'
stop 2
#endif

call enter_exit(sub_name,2)

end subroutine calc_cap_admit_varh

!###################################################################################
!

subroutine Matrix(N_nodes, ha, hv, omega, stiff1, stiff2, RHS1, RHS2)
! This subroutine will create the stiffness matrices and Right hand sides (RHS) for 2 sets of fundamental solutions.
! These matrices is the output of Finite Difference solution of system of ODEs for capillary sheet heights.
!Inputs:
!           1- ha: Non-dimentional sheet height at arterial side of the capillary
!           2- hv: Non-dimentional sheet height at venous side of the capillary
!           3- omega: Dimensionless oscillatory frequency
!           4- N_nodes: The number of nodes in the 1D domain. You can get a better approximation if you increase the number of nodes.
!Outputs:
!           1- Stiff1&2: Stiffness matrices based on different sets of fundamental solutions.
!           2- RHS1&2: Right hand side matrix satisfying the BCs.

! Parameters:
integer, intent(in) :: N_nodes
real(dp), intent(in) :: ha,hv,omega
real(dp), intent(out) :: stiff1(2*N_nodes+2,2*N_nodes+2), stiff2(2*N_nodes+2,2*N_nodes+2) ! Capillary Sheet Stiffness matrices
real(dp), intent(out) :: RHS1(2*N_nodes+2), RHS2(2*N_nodes+2)

! Local Variables
real(dp) :: rep
real(dp) :: dx
real(dp) :: x(N_nodes+1), h(N_nodes+1)
integer :: i ! Matrix row
character(len=60) :: sub_name


sub_name = 'Matrix'
call enter_exit(sub_name,1)


stiff1 = 0 !!! initialise stiff1 matrix
stiff2 = 0 !!! initialise stiff2 matrix
RHS1 = 0 !!! initialise RHS matrix
rep = ha**4 - hv**4
dx = 1.0/(N_nodes - 1.0)  !!! defining delta(x) interval


do i=1,N_nodes+1
    x(i)=(i-1.0)/(N_nodes-1.0)
    h(i)=(ha**4 - rep*x(i))**(0.75)
enddo  !!! x and h^3 array created for x values

    do i = 1,N_nodes+1 !Creating the first stiffness matrix for first fundamental solution + Right hand side
              if (i .eq. N_nodes+1) then
                            stiff1(i,i) = h(i) / (2*dx)
                            stiff1(i, i-2) = -h(i-2) / (2*dx)
                            stiff1(i+N_nodes+1,i+N_nodes+1) = h(i)
                            stiff1(i+N_nodes+1,i+N_nodes-1) = -h(i-2)
                            RHS1(i) = -1.0
              else if (i .eq. N_nodes) then
                            stiff1(i,i) = 1.0
                            stiff1(i+N_nodes+1,i+N_nodes+1) = 1.0
              else
                            stiff1(i,i+2) = h(i+2) / (dx**2)
                            stiff1(i,i+1) = -2.0*h(i+1) / (dx**2)
                            stiff1(i,i) = h(i) / (dx**2)
                            stiff1(i+N_nodes+1,i+N_nodes+3) = h(i+2) / (dx**2)
                            stiff1(i+N_nodes+1,i+N_nodes+2) = -2.0*h(i+1) / (dx**2)
                            stiff1(i+N_nodes+1,i+N_nodes+1) = h(i) / (dx**2)
                            stiff1(i,N_nodes+i+1) = omega
                            stiff1(i+N_nodes+1,i) = -1.0 * omega
              end if
    enddo ! End loop for stiff1

 RHS2 = 0 !!! Reseting right hand side for the next stiffness matrix
 x = 0 !!! Reseting x array for the next stiffness matrix
 h = 0 !!! Reseting h^3 values for the next stiffness matrix

 do i=0,N_nodes
    x(i+1)=(i-1.0)/(N_nodes-1.0)
    h(i+1)=(ha**4 - rep*x(i+1))**(0.75)
enddo  !!! x and h^3 array created for x values


    do i=1,N_nodes+1
            if (i .eq. 1) then
               stiff2(i,i) = -1.0* h(i) / (2*dx)
               stiff2(i,i+2) = h(i+2) / (2*dx)
               stiff2(i+N_nodes+1,i+N_nodes+1) = -1.0 * h(i)
               stiff2(i+N_nodes+1,i+N_nodes+3) = h(i+2)
               RHS2(i) = -1.0
            else if (i .eq. 2) then
               stiff2(i,i) = 1.0
               stiff2(i+N_nodes+1,i+N_nodes+1) = 1.0
            else
               stiff2(i,i-2) = h(i-2) / (dx**2)
               stiff2(i,i-1) = -2.0*h(i-1) / (dx**2)
               stiff2(i,i) = h(i) / (dx**2)
               stiff2(i+N_nodes+1,i+N_nodes-1) = h(i-2) / (dx**2)
               stiff2(i+N_nodes+1,i+N_nodes) = -2.0*h(i-1) / (dx**2)
               stiff2(i+N_nodes+1,i+N_nodes+1) = h(i) / (dx**2)
               stiff2(i,N_nodes+i+1) = omega
               stiff2(i+N_nodes+1,i) = -1.0 * omega
            end if
    enddo
    call enter_exit(sub_name,2)

    end subroutine Matrix

!
!###################################################################################
!

subroutine Mat_to_CC(k,nn,sparsecol,sparserow,sparseval,NonZeros)
! This subroutine forms the sparse representation of matrix K in Compressed Column(CC) format.
! INPUT/OUTPUT(s):
!                - K ====> Stiffness Matrix for Capillary passed from matrix subroutine. (It should be a square matrix)
!                - nn ====> Number of Nodes for the 1D capillary domain. Will be hard coded.
!                - SparseVal ====> An array formed with NonZero values stored in of size(# NonZeros)
!                - SparseRow ====> An array formed with the Row number of each NonZero in matrix of size(# NonZeros)
!                - SparseCol ====> An array formed with values to represent row storage. Includes the index numbers of start of each column (Size of array = MatrixSize + 1)


!Parameters:
real(dp),  intent(in):: k(2*nn+2,2*nn+2)
real(dp),  allocatable, intent(out) :: sparseval(:)
integer, allocatable, intent(out) :: sparsecol(:), sparserow(:)
integer, intent(in) :: nn
integer, intent(out) :: NonZeros
!Local Variables:
integer :: i, j, counter
character(len=60) :: sub_name


sub_name = 'Mat_to_CC'
call enter_exit(sub_name,1)

allocate(sparsecol(2*nn+3)) !allocation of sparsecol
NonZeros = 0 !Initialisation for NonZeros
sparsecol = 1 !SparseCol initialisation

! Counting the number of NonZeros in K matrix
do i = 1,2*nn+2 !going through columns of k
     do j = 1,2*nn+2 !going through rows of k
          if (k(i,j) .NE. 0.0) then
          NonZeros = NonZeros + 1
          end if
     enddo
 enddo
!write(*,*) 'Number of NonZero components in stiffness matrix: ' , NonZeros ! A write statement for checking the number NonZeros

! Allocation of SparseRow and SparseVal (After knowing the number of NonZeros)

 allocate(sparserow(NonZeros))
 allocate(sparseval(NonZeros))

 sparserow = 0 !Initialisation
 sparseval = 0 !Initialisation

 counter=0
 do j = 1,2*nn+2 !going through columns of k
     do i = 1,2*nn+2 !going through rows of k
          if (k(i,j).NE.0) then
            counter = counter + 1
            sparseval(counter)=k(i,j)
            sparserow(counter)=i
          end if
     enddo ! do i
     if (counter.EQ.0.0) then
         sparsecol(j+1) = 0
    else
        sparsecol(j+1) = counter + 1
     end if
 enddo ! do j
call enter_exit(sub_name,2)

end subroutine Mat_to_CC

end module capillaryflow
