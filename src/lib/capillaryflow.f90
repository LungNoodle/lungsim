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
  use other_consts, only: TOLERANCE
  implicit none
  

  !Module parameters
  
  !Module types

  !Module variables

  !Interfaces
  private
  public cap_flow_ladder
  
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

    use diagnostics, only: enter_exit
    use arrays, only:dp,capillary_bf_parameters,num_units
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
    use diagnostics, only: enter_exit
    use solve, only: pmgmres_ilu_cr
    use arrays, only:dp,capillary_bf_parameters

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
    use diagnostics, only: enter_exit
    use arrays, only:dp,capillary_bf_parameters

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
    use diagnostics, only: enter_exit
    use arrays, only:dp,capillary_bf_parameters
    use other_consts, only:PI

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
    use arrays, only: dp
    use diagnostics, only: enter_exit
    use arrays, only:dp,capillary_bf_parameters

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
    use arrays, only: dp
    use diagnostics, only: enter_exit
    use arrays, only:dp,capillary_bf_parameters

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
!...  Stepping down linearly with each generation from Fung's
!...  estimate at 45% hematocrit (4e-3Pa.s) to his estimate at 30% hematocrit
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

    call enter_exit(sub_name,2)
end subroutine cap_specific_parameters

end module capillaryflow

