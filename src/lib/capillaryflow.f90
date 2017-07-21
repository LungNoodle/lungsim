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
  use other_consts
  implicit none
  
  !Module parameters
  
  !Module types

  !Module variables

  !Interfaces
  private
  public cap_flow_ladder,cap_flow_sheet,evaluate_ladder,read_capillary_params
  
contains
!
!###################################################################################
!

  subroutine cap_flow_ladder(ne,HEIGHT,LPM_FUNC,Pin,Pout,Ppl,&
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

    implicit none
    integer :: ne
    real(dp) :: HEIGHT(3),LPM_FUNC,Pin,Pout,Ppl,Q01,R_in,R_out,x,y,z
    logical :: OUTPUT_PERFUSION
    
    !    Local variables
    integer :: MatrixSize,NonZeros,submatrixsize,i,N
    real(dp) :: area,Pressure(4*num_symm_gen-4),Q01_mthrees,sheet_number
    character(len=60) :: sub_name
    
    sub_name = 'cap_flow_ladder'
    call enter_exit(sub_name,1)
    
    !CALL ASSERT(num_symm_gen.LT.13,'>>max symmetrics [13] exceeded!',
    ! &   ERROR,*9999)

!     Number of non-zero entries in solution matrix. 
      NonZeros=3
      DO i=2,num_symm_gen
         NonZeros=NonZeros+4*i+10
      ENDDO
!     The size of the solution matrix (number of unknown pressures and flows)
      MatrixSize=5*num_symm_gen-3
!     The number of unknown pressures
      submatrixsize=4*num_symm_gen-4
  
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
      DO i=1,num_symm_gen
         sheet_number=sheet_number+2.d0**i
      ENDDO
      area=total_cap_area/sheet_number !m^2
!...  Initial guess for pressure distribution lets say all arterial pressures are the same 
!...  and all the venous pressures are the same solution appears independent of this.
      DO i=1,num_symm_gen-1
          Pressure(4*i-3)=1000 ! Pa
          Pressure(4*i-2)=1000
          Pressure(4*i-1)=100
          Pressure(4*i)=100
      ENDDO
     
!  ---CALL THE FUNCTIONS THAT CALCULATE THE FLOW ACROSS THE LADDER FOR A GIVEN PRESSURE DROP--
      call evaluate_ladder(ne,NonZeros,MatrixSize,submatrixsize,&
      area,HEIGHT,Pin,Pout,Ppl,Pressure,Q01_mthrees,x,y,z,&
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
!###################################################################################
!
  subroutine cap_flow_sheet(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT,&
     zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited)
!*cap_flow_sheet:* Calculates resistance across a capillary sheet.
    use diagnostics, only: enter_exit
    use other_consts
    implicit none
    integer :: ne,zone
    real(dp) :: SHEET_RES,Q_c,Hart,Hven,Pin_sheet,Pout_sheet,Ppl,area
    
    !     Local Variables
    real(dp) :: area_new,C,Hmax_art,Hmax_ven,Hub,L_new,P_a,Palv_pa,P_ave,&
         P_v,RBC_tt,recruited,waterfall_scale
    character(len=60) :: sub_name
    
    sub_name = 'cap_flow_sheet'
    call enter_exit(sub_name,1)
    
    !...  area of sheet is an input. Now area and pathlength need to be scaled. 
    area_new=area*area_scale
    L_new=L_c*L_scale
    !... Maximum sheet heights
    Hmax_art=H0+alpha_c*Pub_c 
    Hmax_ven=Hmax_art
    !...  Blood pressure in and out of capillary sheet
    P_a=Pin_sheet           
    P_v=Pout_sheet
    Palv_pa=Palv*98.06d0
    !...   CAPILLARY RECRUITMENT
    !...   CALCULATE AVERAGE CAPILLARY BLOOD PRESSURE
    P_ave=(P_a+P_v)/2.d0
    !      calculate recruitment proportion
    recruited=1.d0-F_rec*EXP(-(P_ave**2.d0)/(sigma_rec**2.d0))
    !      multiply flows and transit times
    area_new=area_new*recruited
    
    IF((P_v-Palv_pa).LE.Plb_c)THEN
       ! ### Zone 2 (waterfall) scaling 
       waterfall_scale=1-F_sheet+F_sheet* &
          EXP(-(P_v-Palv_pa)**2.d0/(2.d0*sigma_cap**2.d0))
       area_new=area_new*waterfall_scale
    ENDIF
    !  -----CALCULATING SHEET FLOW----
    !...  Sheet flow model constant, C
    C=area_new/(4.d0*mu_c*K_cap*F_cap*L_new**2*alpha_c)
    !###  Determine what zone we are in and calculate flow and TT
    IF((P_a-Palv_pa).LT.Plb_c)THEN 
       !...   ZONE 1:
       !...   Arteriole and venous pressure both less than alveolar pressure - the
       !...   sheet is collapsed
       zone=1
       Hart=0.d0
       Hven=0.d0            
       Q_c=0.d0
       RBC_tt=0.d0 !TEMPORARY
    ELSEIF((P_a-Palv_pa).LE.Pub_c.AND.(P_v-Palv_pa).LE.Plb_c)THEN
       !...    ZONE 2:
       zone=2
       Hart=H0+alpha_c*(P_a-Palv_pa)
       Hven=H0
       Q_c=C*(Hart**4.d0-H0**4.d0)
       RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c) &
            /(Hart**3.d0-H0**3.d0)
    ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).LE.Plb_c)THEN
       !...       ZONE 2:
       zone=2
       Hart=Hmax_art
       Hven=H0
       Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-Palv_pa-Pub_c) &
            +(Hmax_art**4.d0-H0**4.d0)/(4.d0*alpha_c))
       RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c)/ &
            ((3.d0*alpha_c*Hmax_art**2.d0* &
            (P_a-Palv_pa-Pub_c)+(Hmax_art**3.d0-H0**3.d0))) 
    ELSEIF((P_a-Palv_pa).LE.Pub_c.AND.(P_v-Palv_pa).LE.Pub_c)THEN
       !...       ZONE3 or 4: The boundary between zone 3 and 4 is not clearcut.
       zone=3 !tmp!!should = 3
       Hart=H0+alpha_c*(P_a-Palv_pa)
       Hven=H0+alpha_c*(P_v-Palv_pa)
       Q_c=C*(Hart**4.d0-Hven**4.d0)
       RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c) &
            /(Hart**3.d0-Hven**3.d0)
    ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).LE.Pub_c)THEN
       zone=3
       Hart=Hmax_art
       Hven=H0+alpha_c*(P_v-Palv_pa)
       Hub=H0+alpha_c*Pub_c
       Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-Palv_pa-Pub_c)+ &
            (Hub**4.d0-Hven**4.d0)/(4*alpha_c))
       RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2*alpha_c) &
            /(3.d0*alpha_c*Hmax_art**2.d0*(P_a-Palv_pa-Pub_c)+ &
            (Hmax_art**3.d0-Hven**3.d0))
    ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).GT.Pub_c)THEN
       zone=4 !!!tmp should = 3 
       Hart=Hmax_art
       Hven=Hmax_ven
       Q_c=4.d0*C*alpha_c*Hart**3.d0*(P_a-P_v)
       RBC_tt=(mu_c*F_cap*K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
    ELSE
       !.... SOMETHING HAS GONE TERRIBLY WRONG AS THIS SHOULD BE ALL THE OPTIONS!!
       Q_c=1.d-8
       RBC_tt=0.d0
       write(*,*) 'Error, incorrect option, something is wrong!'
    ENDIF

    !... RBC transit time (1.4 times faster than blood)
    RBC_tt=RBC_tt/1.4d0
    IF(alpha_c.eq.0.d0)THEN
       C=area_new/(4.d0*mu_c*K_cap*F_cap*L_new**2)
       Hart=Hmax_art
       Hven=Hmax_ven
       Q_c=4.d0*C*Hart**3.d0*(P_a-P_v)
       RBC_tt=(mu_c*F_cap*K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
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
!######################################################################################################
!

  subroutine evaluate_ladder(ne,NonZeros,MatrixSize,submatrixsize, &
      area,HEIGHT,Pin,Pout,Ppl,Pressure,Q01_mthrees,x,y,z, &
      OUTPUT_PERFUSION)
!*evaluate_ladder:* Sets up and solves matrix equations for ladder model.
    use diagnostics, only: enter_exit
    use solve, only: BICGSTAB_LinSolv
    implicit none
    integer :: ne,NonZeros,MatrixSize,submatrixsize
    real(dp) :: area,HEIGHT(3),Pin,Pout,Ppl,Pressure(submatrixsize),&
      Q01_mthrees,x,y,z
    logical :: OUTPUT_PERFUSION

    ! Local variables
    integer :: i,iter,j,gen,SparseCol(NonZeros),SparseRow(MatrixSize+1),&
      zone,nj,num_sheet
    real(dp) :: area_new,ErrorEstimate,Hart,Hven,Pin_SHEET,Pout_SHEET,&
         Q_c,Qtot,Qgen,Q_sheet(num_symm_gen),RBC_TT,RHS(MatrixSize),&
         Rtot,SHEET_RES,Solution(MatrixSize),SolutionLast(MatrixSize),&
         SparseVal(Nonzeros),TOTAL_CAP_VOL,TOTAL_SHEET_H,TOTAL_SHEET_SA,&
         P_inlet,P_outlet,R_upstream,R_downstream,GRAVITY,recruited,TT_TOTAL
    character(len=60) :: sub_name

    sub_name = 'evaluate_ladder'
    call enter_exit(sub_name,1)

!###  INITIAL SOLUTION GUESS 
      DO i=1,submatrixsize
      solution(i)=Pressure(i)
      ENDDO
      DO i=submatrixsize+1,matrixSize-1
       solution(i)=Q01_mthrees/2**num_symm_gen
      ENDDO
      Solution(Matrixsize)=Q01_mthrees
!###  INITIALISE SOLUTIONLAST
      DO j=1,MatrixSize
        SolutionLast(j)=Solution(j)
      ENDDO
     
!### INPUT TO THE LADDER MODEL THAT IS INDEPENDENT OF ITERATION
      CALL LADDERSOL_MATRIX(NonZeros,MatrixSize,submatrixsize,&
         SparseCol,SparseRow,SparseVal,RHS,Pin,Pout)

!### ITERATIVE LOOP
      iter=0 
      ErrorEstimate=1.d10
      DO WHILE(ErrorEstimate.GT.1.0d-9.AND.iter.LT.100)
        iter=iter+1
!...  CALCULATE RESISTANCE GIVEN CURRENT PRESSURE AND FLOW - THEN UPDATE 
!...  SparseVal- THese are the only elements of the solution matrix that need
!.... iteratively updating
        CALL POPULATE_MATRIX_LADDER(ne,NonZeros,MatrixSize,&
         submatrixsize,SparseCol,SparseRow,area,Pin,Pout,Ppl,Pressure,&
         Q01_mthrees,Q_sheet,RHS,SparseVal)
 
!        CALL MINI_LINEAR_SYSTEM(MatrixSize,NonZeros,SparseCol,SparseRow,&
!        Solution,SparseVal,RHS)

        call BICGSTAB_LinSolv(MatrixSize,NonZeros,RHS,Solution,SparseCol,&
             SparseRow,SparseVal,,1.d-9,1000) !NB/ 1.d-9 = convergence tol, 1000 = max solver iterations

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
      Do i=1,num_symm_gen
        Qtot=Qtot+Q_sheet(i)*2.d0**i
      ENDDO
      Rtot=(Pin-Pout)/Q01_mthrees
       
       IF(OUTPUT_PERFUSION)THEN
!###  GET SOLUTIONS TO WRITE TO FILE 
        TOTAL_CAP_VOL=0.d0
        TOTAL_SHEET_SA=0.d0
        TT_TOTAL=0.d0
        TOTAL_SHEET_H=0.d0
        num_sheet=0
!... Ladder output
        DO i=1,num_symm_gen-1
           gen=i
           Pin_sheet=Pressure(4*i-3)
           Pout_sheet=Pressure(4*i-1)
!...     calulate resistance and transit time through a single capillary
           CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
           zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited)
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

        gen=num_symm_gen 
        Pin_sheet=Pressure(4*num_symm_gen-6) !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*num_symm_gen-4) !%pressure out of final capillary sheets
         CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
              zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited)
         num_sheet=num_sheet+2**gen
           Qgen=Q_C*2.d0**num_symm_gen
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

    call enter_exit(sub_name,2)

  end subroutine evaluate_ladder

!
!######################################################################################################
!

subroutine laddersol_matrix(NonZeros,MatrixSize,submatrixsize,&
      SparseCol,SparseRow,SparseVal,RHS,Pin,Pout)
!*laddersol_matrix:*Sets up ladder matrix entries that are independent of iteration.
    use diagnostics, only: enter_exit

    implicit none
    integer :: NonZeros,MatrixSize,submatrixsize,SparseRow(MatrixSize+1),&
       SparseCol(Nonzeros)
    real(dp) :: SparseVal(Nonzeros),RHS(MatrixSize),Pin,Pout
    !    Local variables
    integer :: count1,count,i,j
    character(len=60) :: sub_name
    
    sub_name = 'cap_flow_ladder'
    call enter_exit(sub_name,1)

!  Define the vector Ap (for UMFPACK) this vector starts at zero and then each entry is the cumulative total of non-zero entries as you step through columns from L to right
      SparseRow(1)=1
      SparseRow(2)=3
      SparseRow(3)=7
      SparseRow(4)=9
      SparseRow(5)=13
      DO i=2,num_symm_gen-1
         SparseRow(6+4*(i-2))=SparseRow(6+4*(i-2)-1)+2+i
         SparseRow(7+4*(i-2))=SparseRow(6+4*(i-2))+3+i
         SparseRow(8+4*(i-2))=SparseRow(7+4*(i-2))+2+i
         SparseRow(9+4*(i-2))=SparseRow(8+4*(i-2))+3+i
      ENDDO
      DO i=1,num_symm_gen
        SparseRow(4*num_symm_gen-3+i)=SparseRow(submatrixsize+i)+3
      ENDDO
      SparseRow(MatrixSize+1)=SparseRow(MatrixSize)+num_symm_gen+1
      
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
      DO i=1,num_symm_gen
         SparseCol(NonZeros-i)=MatrixSize-i
         SparseVal(NonZeros-i)=-2.d0**(num_symm_gen+1-i)
         SparseCol(NonZeros-num_symm_gen-1-3*(i-1))=MatrixSize-i
      ENDDO
      !...  CAPILLARY CONNECTIONS
      !...  Prior to final generation..      
      DO i=1,num_symm_gen-1
        SparseCol(NonZeros-4*num_symm_gen+3*(i-1))=1+4*(i-1)
        SparseVal(NonZeros-4*num_symm_gen+3*(i-1))=1.d0
        SparseCol(NonZeros-4*num_symm_gen+3*(i-1)+1)=3+4*(i-1)
        SparseVal(NonZeros-4*num_symm_gen+3*(i-1)+1)=-1.d0
      ENDDO
      !...  FINAL GENERATION
      SparseCol(NonZeros-num_symm_gen-2)=MatrixSize-num_symm_gen-1
      SparseVal(NonZeros-num_symm_gen-2)=-1.d0
      SparseCol(NonZeros-num_symm_gen-3)=MatrixSize-num_symm_gen-3
      SparseVal(NonZeros-num_symm_gen-3)=1.d0
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
      DO i=2,num_symm_gen-1
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
subroutine populate_matrix_ladder(ne,NonZeros,MatrixSize,&
      submatrixsize,SparseCol,SparseRow,area,Pin,Pout,Ppl,Pressure,&
      Q01_mthrees,Q_sheet,RHS,SparseVal)
!*populate_matrix_ladder:*Sets up ladder matrix entries that are NOT independent of iteration.
    use diagnostics, only: enter_exit
    use other_consts

    implicit none
    integer :: ne,NonZeros,MatrixSize,submatrixsize,SparseCol(NonZeros),&
           SparseRow(MatrixSize+1)
    real(dp) :: area,Pin,Pout,Ppl,Pressure(submatrixsize),Q01_mthrees,&
            Q_sheet(num_symm_gen),RHS(MatrixSize),SparseVal(Nonzeros)
    !    Local variables
    integer :: gen,count,j,zone
    real(dp) :: Q,Palv_Pa,P_exta,P_extv,radupdate,R_art1,R_art2,&
        R_ven1,R_ven2,SHEET_RES,Q_c,R_sheet(num_symm_gen),Q_gen,Hart,&
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
      DO gen=1,num_symm_gen-1
!...    FIRST HALF OF ARTERIOLE
!...    Update radius of arteriole based on inlet pressure 
         IF(rad_a(gen).LT.100.d-6) THEN 
           P_exta=Palv*98.06d0 ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_exta=-Ppl
         ENDIF
         IF ((Pressure(4*gen-3)-P_exta).LE.Pub_a_v)THEN
!            radupdate=rad_a(gen)*(1+(Pressure(4*gen-3)-P_exta)*alpha_a)
            radupdate=rad_a(gen)+alpha_a*(Pressure(4*gen-3)-P_exta)*&
            (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-3)-P_exta)/(2.d0-2.d0*num_symm_gen)
         ELSE
!            radupdate=rad_a(gen)*(1+Pub_a_v*alpha_a)
            radupdate=rad_a(gen)+alpha_a*Pub_a_v*(gen-num_symm_gen)&
            /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v&
            /(2.d0-2.d0*num_symm_gen)
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
           IF (Pressure(4*gen-1)-Palv_pa.LE.Pub_a_v)THEN
           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-Palv_pa)*alpha_v)
           ELSE
              radupdate=rad_v(gen)*(1+Pub_a_v*alpha_v)
           ENDIF
         IF(rad_v(gen).LT.100.d-6) THEN 
           P_extv=Palv*98.06d0 ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_extv=-Ppl
         ENDIF
         IF((Pressure(4*gen-1)-P_extv).LE.Pub_a_v)THEN
!           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-P_extv)*alpha_v)
            radupdate=rad_v(gen)+alpha_v*(Pressure(4*gen-1)-P_extv)*&
            (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-1)-P_extv)/(2.d0-2.d0*num_symm_gen)
         ELSE
!            rad_update=rad_v(gen)*(1+Pub_a_v*alpha_v)
            radupdate=rad_v(gen)+alpha_v*Pub_a_v*(gen-num_symm_gen)&
            /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v&
            /(2.d0-2.d0*num_symm_gen)
         ENDIF
         volume_vessels=volume_vessels+&
          (2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0)

!...   Calculate Poiseuille resistance in first half of venule
           R_ven1=(8.d0*mu_app(gen)*L_v(gen)/2.d0)/(pi*radupdate**4.d0);

!...   CAPILLARY ELEMENT (arteriole + venule + capillary)  
!...    pressure into the capillaries
           Pin_sheet=Pressure(4*gen-3)
           Pout_sheet=Pressure(4*gen-1)
!...     calulate resistance and transit time through a single capillary
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited)
            Q_sheet(gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(gen)*2**gen
            R_sheet(gen)=SHEET_RES
!...    Update soln matrix (delta p=QR across a cap)
           SparseVal(NonZeros-num_symm_gen-1-3*(num_symm_gen-gen))=&
            -SHEET_RES
 
     !...   SECOND HALF OF ARTERIOLE
!...    Update radius of arteriole based on inlet pressure NB: May need to
!....   update this later to give an average of inlet and outlet radii
      IF (Pressure(4*gen-2)-P_exta.LE.Pub_a_v)THEN
!            radupdate=rad_a(gen)*(1+(Pressure(4*gen-2)-P_exta)*alpha_a)
            radupdate=rad_a(gen)+alpha_a*(Pressure(4*gen-2)-P_exta)*&
            (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen-2)-P_exta)/(2.d0-2.d0*num_symm_gen)
         ELSE
!            radupdate=rad_a(gen)*(1+Pub_a_v*alpha_a)
            radupdate=rad_a(gen)+alpha_a*Pub_a_v*(gen-num_symm_gen)&
            /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v&
            /(2.d0-2.d0*num_symm_gen)
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
      IF (Pressure(4*gen)-P_extv.LE.Pub_a_v)THEN
!           radupdate=rad_v(gen)*(1+(Pressure(4*gen)-P_extv)*alpha_v)
            radupdate=rad_v(gen)+alpha_v*(Pressure(4*gen)-P_extv)*&
            (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*&
            (Pressure(4*gen)-P_extv)/(2.d0-2.d0*num_symm_gen)
         ELSE
!            rad_update=rad_v(gen)*(1+Pub_a_v*alpha_v)
            radupdate=rad_v(gen)+alpha_v*Pub_a_v*(gen-num_symm_gen)&
            /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v&
            /(2.d0-2.d0*num_symm_gen)
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
        Pin_sheet=Pressure(4*num_symm_gen-6); !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*num_symm_gen-4); !%pressure out of final capillary sheets
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, &
        zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited)
        SparseVal(NonZeros-num_symm_gen-1)=-SHEET_RES
            Q_sheet(num_symm_gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(num_symm_gen)*2**gen
            R_sheet(num_symm_gen)=SHEET_RES

    call enter_exit(sub_name,2)

end subroutine populate_matrix_ladder

!
!######################################################################################################
!

subroutine read_capillary_params
   use other_consts
   use arrays, only : num_units
   implicit none

  ! Input related variables
  character(len=100) :: buffer, label
  integer :: pos
  integer, parameter :: fh = 15
  integer :: ios = 0
  integer :: line = 0

  ! Control file variables
!  real :: density,viscosity,Go,beta,Go_vein,beta_vein
!  integer, dimension(5) :: vector

  open(fh, file='Parameters/capillary.txt')

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('Sheet_height')
           read(buffer, *, iostat=ios) sheet_h0
           !!print *, 'Read sheet height: ', sheet_h0
        case ('K_factor')
           read(buffer, *, iostat=ios) K_cap
           !print *, 'Read K friction factor: ', K_cap
        case ('F_factor')
           read(buffer, *, iostat=ios) F_cap
           !print *, 'Read F friction factor: ', F_cap
        case ('Zone2_F_constant')
           read(buffer, *, iostat=ios) F_sheet
           !print *, 'Read F_sheet: ', F_sheet
        case ('Zone2_sigma_constant')
           read(buffer, *, iostat=ios) sigma_cap
           !print *, 'Read sigma_cap: ', sigma_cap
        case ('Apparent_viscosity')
           read(buffer, *, iostat=ios) mu_c
           !print *, 'Read mu_c: ', mu_c
         case ('Recruitment_parameters_F')
           read(buffer, *, iostat=ios) F_rec
           !print *, 'Read F_rec: ', F_rec
         case ('Recruitment_parameters_sigma')
           read(buffer, *, iostat=ios) sigma_rec
           !print *, 'Read sigma_rec: ', sigma_rec
         case ('Path_length_art_ven')
           read(buffer, *, iostat=ios) L_c
           !print *, 'Read L_c: ', L_c
         case ('Lower_pressure_sheet')
           read(buffer, *, iostat=ios) Plb_c
           !print *, 'Read Plb_c: ', Plb_c
         case ('Upper_pressure_sheet')
           read(buffer, *, iostat=ios) Pub_c
           !print *, 'Read Pub_c: ', Pub_c
         case ('Upper_pressure_acinar')
           read(buffer, *, iostat=ios) Pub_a_v
           !print *, 'Read Pub_a_v: ', Pub_a_v
         case ('Capillary_surface_area')
           read(buffer, *, iostat=ios) total_cap_area
           !print *, 'Read total_cap_area: ', total_cap_area
           IF(num_units.GT.0)THEN
              total_cap_area=total_cap_area/num_units
           ELSE       
              write(*,*) 'WARNING: NTERMINAL is zero - assuming 32000 acini'
              write(*,*) 'remember to do fem evaluate order!'
              total_cap_area=total_cap_area/32000d0
           ENDIF
                  
         case ('Terminal_arteriole_length')
           read(buffer, *, iostat=ios) L_art_terminal
           !print *, 'Read Terminal_arteriole_length:', L_art_terminal
         case ('Terminal_venule_length')
           read(buffer, *, iostat=ios) L_vein_terminal
           !print *, 'Read Terminal_venule_length: ', L_vein_terminal
         case ('Terminal_arteriole_radius')
           read(buffer, *, iostat=ios) R_art_terminal
           !print *, 'Read Terminal_arteriole_radius: ', R_art_terminal
         case ('Terminal_venule_radius')
           read(buffer, *, iostat=ios) R_vein_terminal
           !print *, 'Read Terminal_venule_radius: ', R_vein_terminal
         case ('Initial_acinar_model_resistance')
           read(buffer, *, iostat=ios) INITIAL_LPM
           !print *, 'Read Initial_acinar_model_resistance: ', INITIAL_LPM
         case ('Upper_pressure_bound_for_extra-acinar_vessels')
           read(buffer, *, iostat=ios) Ptm_max
           !print *, 'Read Upper_pressure_bound_for_extra-acinar_vessels: ', Ptm_max
         case ('Number_of_symmetric_generations')
           read(buffer, *, iostat=ios) num_symm_gen
           !print *, 'Read Number_of_symmetric_generations: ', num_symm_gen
         case default
           print *, 'Skipping invalid label at line', line
        end select
     end if
  end do
  
end subroutine read_capillary_params


!
!######################################################################################################
!
end module capillaryflow

