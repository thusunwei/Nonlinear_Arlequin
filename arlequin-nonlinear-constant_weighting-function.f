C###################################################################################################
C  USER DEFINED ELEMENT:2D 4NODES ISOPARAMETER LINEAR-ELASTIC MATERIA 
C  THIS CODE IS PROGRAMMED IN ONE PURPOSED: 
C  To implement complement the combination of UEL and UMAT     
C--------------------------------------------------------------------------------------------------C
C  MOST GENERAL  2D STRUCTURAL ELEMENT HAS TWO MATERIAL PROPERTIES: E  THICK V 
C  STATE VARIABLES ARE STRESS AND STRAIN IN EACH POI        
C  NSVARS:4(INTEGRATION POINT)X(7)=28,THEN SVARS   
C  NDOFEL:THE DEGREE OF FREEDOM= 4*2    
C###################################################################################################
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(*),ENERGY(7),PROPS(3),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(4),
     5     JPROPS(*)
C
C
C***** DEFINE  PARAMETERS AND M     
    
    
      PARAMETER (HALF=0.5D0,ZERO=0.D0,ONE=1.D0,TWO=2.D0,
     1  LLL=1)   !need to be corrected
      PARAMETER (THREE=3.0D0)    
      DIMENSION Gauss(4,2)                                        !高斯积分点坐标
      DIMENSION BN(2,4),BNV(2,4)                                  !N矩阵偏导
      DIMENSION P(2,8),PV(2,8)                                    !N矩阵
      DIMENSION BJ(2,2),BL(2,4),BR(2,2)                           !雅克比矩阵和过渡矩阵
      DIMENSION BJV(2,2),BLV(2,4),BRV(2,2)                        !虚拟单元雅克比矩阵和过渡矩阵                    
      DIMENSION B(3,8),D(3,3),S(3,8)                              !B阵、D阵、S阵
      DIMENSION BV(3,8),DV(3,3),SV(3,8)                           !虚拟单元B阵、D阵、S阵
      DIMENSION SSTRAIN(3),SSTRESS(3)                             !应力应变矩阵
      DIMENSION CVJ(8,8),PVNJ(8,8),S1(8,8)                        !应力应变矩阵
      DIMENSION BVBJ(8,8),CVJTRANSA(8,8),AMA(8,8)                 !应力应变矩阵
      DIMENSION BVTRANSA(8,3),PVTRANSA(8,2),BTRANSA(8,3)          !应力应变矩阵
      DIMENSION BODYFORCE(8,1),CDF_D(4,4)  
      DIMENSION CDF_STRAIN(4),CDF_DSTRAIN(4),CDF_STRESS(4)          
      DIMENSION CDF_EELAS(4),CDF_EPLAS(4),UMAT_SVARS(11) 
      DIMENSION CDF_TSTRESS(3)
      PARAMETER (nelem=1000000,ninpt=4,nsvint=19)    
      DIMENSION temporary_array(nsvint)
      common /kuser/user_variable(nelem,nsvint,ninpt)
      
      if (JTYPE.EQ.1) THEN    
     
      Gauss(1,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(1,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(2,1)=SQRT(ONE/THREE)
      Gauss(2,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(3,1)=SQRT(ONE/THREE)
      Gauss(3,2)=SQRT(ONE/THREE)
      Gauss(4,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(4,2)=SQRT(ONE/THREE)*1.0D0
C      
C------------------------------------
C    INITAIL  RHS、AMATRX
C    NDOFEL：节点的自由度数 
C    NRHS：  Number of load vectors.
C    AMATRX: 刚度矩阵
C    Gauss： 高斯积分
C----------------------------------
      DO K1 =1,NDOFEL
           DO KRHS = 1, NRHS
             RHS(K1,KRHS) = ZERO
           EN DDO
           DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
           ENDDO
      ENDDO
      

      DO Kintk=1,NNODE
    
           xi=Gauss(Kintk,1)
           eta=Gauss(Kintk,2)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Make B Matrix  in Natural Coordinate System                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****BN:Shape Functions     
           DO I=1,2
             DO J=1,NNODE
               BN(I,J)=ZERO
             ENDDO
           ENDDO           
C*****N为形函数矩阵
            DO I=1,2
             DO J=1,8
               P(I,J)=ZERO
             ENDDO
           ENDDO
      P(1,1)=0.25*(1.0*ONE-ETA)*(1.0*ONE-XI)
      P(1,3)=0.25*(ONE-ETA)*(ONE+XI)
      P(1,5)=0.25*(ONE+ETA)*(ONE+XI)
      P(1,7)=0.25*(1.0*ONE+ETA)*(ONE-XI)
      P(2,2)=0.25*(1.0*ONE-ETA)*(1.0*ONE-XI)
      P(2,4)=0.25*(ONE-ETA)*(ONE+XI)
      P(2,6)=0.25*(ONE+ETA)*(ONE+XI)
      P(2,8)=0.25*(1.0*ONE+ETA)*(ONE-XI)      
C*****BJ:Jacobi Matrix
           DO I=1,2
             DO K=1,2
                BJ(I,K)=0.0D0
             ENDDO
           ENDDO
           
          CALL F_YAKEBI(BJ,BN,COORDS,xi,eta,4)      !计算雅可比矩阵
          
C*****BR Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BR(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BH,BJ)
           
           BR(1,1)=BJ(2,2)/BH
           BR(1,2)=-1*BJ(1,2)/BH
           BR(2,1)=-1*BJ(2,1)/BH
           BR(2,2)=BJ(1,1)/BH
           
C           
C*****BL: Transition Matrix
           DO I=1,2
             DO K=1,NNODE
               BL(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BL,BR,2,2,BN,2,4)  !BRxBN
C*****Make B the Matrix         
           DO I=1,3
             DO K=1,8
               B(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,NNODE
             B(1,K*2-1)=BL(1,K)
             B(2,K*2)=BL(2,K)
             B(3,K*2-1)=BL(2,K)
             B(3,K*2)=BL(1,K)
           ENDDO           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE STRESSES AND STRAINS                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****INITIALIZE CDF_STRAIN,CDF_DSTRAIN, CDF_STRESS
           DO I=1,4
              CDF_STRAIN(I)=ZERO
              CDF_DSTRAIN(I)=ZERO
              CDF_STRESS(I)=ZERO
           ENDDO           
          
           DO K1=1,3
             DO K2=1,3
              D(K1,K2)=ZERO
             END DO
           END DO  	           
C
C*****由位移和位移增量求得应变和应变增量，同时把3X1向量，拓展到4X1
C     
           DO I=1,3
              DO J=1,8
                CDF_STRAIN(I)=CDF_STRAIN(I)+B(I,J)*U(J)
                CDF_DSTRAIN(I)=CDF_DSTRAIN(I)+B(I,J)*DU(J,1)
              ENDDO
           ENDDO 
C           
           CDF_STRAIN(4)=CDF_STRAIN(3)  
           CDF_STRAIN(3)=ZERO
           CDF_DSTRAIN(4)=CDF_DSTRAIN(3)  
           CDF_DSTRAIN(3)=ZERO
C
C*****从上一载荷步载入塑性应变、应力以及等效塑性应变的值
          
         DO K1=1,11
               UMAT_SVARS(K1)=SVARS(11*Kintk+K1-11)
          END DO     
          
           	
C*****调用子函数F_UMAT材料矩阵以及输出应变、应力      
           CALL F_UMAT_MISES(CDF_STRAIN,CDF_DSTRAIN,NPROPS,PROPS,
     1            CDF_D,UMAT_SVARS)        
C
C***** MAKE S-MATRIX  
C  
C*****提取平面应变D材料矩阵，构建平面应变S阵                     
            DO I=1,3
                DO K=1,8
                   S(I,K)=0.0D0
              ENDDO
           ENDDO
           
           DO K1=1,2
             DO K2=1,2
              D(K1,K2)=CDF_D(K1,K2)
             END DO
              D(3,K1)=CDF_D(4,K1) 	
           END DO 
              D(3,3)=CDF_D(4,4)      
 
              DO I=1,3
                DO K=1,8
                   DO J=1,3
                      S(I,K)=S(I,K)+D(I,J)*B(J,K)
                   ENDDO
              ENDDO
           ENDDO       
C*****提取应力构建平面应变应力向量（3X1）             
           DO I=1,2
             CDF_TSTRESS(I)=UMAT_SVARS(I+4)
           ENDDO 
             CDF_TSTRESS(3)=UMAT_SVARS(8)       
             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE body force                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
         DO K1=1,8
		    BODYFORCE(K1,1)=0
         ENDDO
        	    BODYFORCE(2,1)=-P(1,1)*0*10*BH	
		    BODYFORCE(4,1)=-P(1,3)*0*10*BH
		    BODYFORCE(6,1)=-P(1,5)*0*10*BH
		    BODYFORCE(8,1)=-P(1,7)*0*10*BH	  
              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                     ASSEMBLE RHS AND LHS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ASSEMBLE RHS AND LHS

           CALL F_TRANS(BTRANSA,B,3,8)
            
          DO K1=1,8
          RHS(K1,1)=RHS(K1,1)+BODYFORCE(K1,1)
          END DO

           DO K1=1,8
			    DO KK=1,3
                 RHS(K1,1)=RHS(K1,1)-BTRANSA(K1,KK)*CDF_TSTRESS(KK)*BH
                 ENDDO
          ENDDO
          
          DO K1=1,8
              DO K2=1,8
                 DO KK=1,3
              AMATRX(K1,K2)=AMATRX(K1,K2)+BTRANSA(K1,KK)*S(KK,K2)*BH
                 ENDDO
              ENDDO
           ENDDO
C                      
C***** STORE STRAINS IN STATE VARIABLE ARRAY
C
           DO K1=1,11
               SVARS(11*Kintk+K1-11)=UMAT_SVARS(K1)
           END DO     
           
           DO I=1,11
              temporary_array(I)=UMAT_SVARS(I)  
          ENDDO
          
         if(kintk.eq.1)then
           do i=1,11
           user_variable(jelem,I,Kintk)=temporary_array(I)
           enddo
           endif
           
           if(kintk.eq.2)then
            do i=1,11
           user_variable(jelem,I,Kintk)=temporary_array(I)
           enddo
           endif
    
           if(kintk.eq.3)then
             do i=1,11
           user_variable(jelem,I,Kintk+1)=temporary_array(I)
           enddo
           endif
           
           if(kintk.eq.4)then
             do i=1,11
           user_variable(jelem,I,Kintk-1)=temporary_array(I)
           enddo
           endif                 

      ENDDO
     
C*****guass point end 
      
      ENDIF
      
      
C      
C------------------------------------
C    粗糙单元与虚拟单元耦合
C----------------------------------      
       
      if (JTYPE.EQ.2) THEN    
      
      Gauss(1,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(1,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(2,1)=SQRT(ONE/THREE)
      Gauss(2,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(3,1)=SQRT(ONE/THREE)
      Gauss(3,2)=SQRT(ONE/THREE)
      Gauss(4,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(4,2)=SQRT(ONE/THREE)*1.0D0
C      
C------------------------------------
C    INITAIL  RHS、AMATRX
C    NDOFEL：节点的自由度数 
C    NRHS：  Number of load vectors.
C    AMATRX: 刚度矩阵
C    Gauss： 高斯积分
C----------------------------------
      DO K1 = 1, NDOFEL
           DO KRHS = 1, NRHS
             RHS(K1,KRHS) = ZERO
           ENDDO
           DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
           ENDDO
      ENDDO
     
       DO Kintk=1,4
           xi=Gauss(Kintk,1)
           eta=Gauss(Kintk,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Make B Matrix  in Natural Coordinate System                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****BN:Shape Functions     
           DO I=1,2
             DO J=1,4
               BN(I,J)=ZERO
             ENDDO
           ENDDO           
C*****N为形函数矩阵
            DO I=1,2
             DO J=1,8
               P(I,J)=ZERO
             ENDDO
           ENDDO
      P(1,1)=0.25D0*(1.0D0*ONE-ETA)*(1.0D0*ONE-XI)
      P(1,3)=0.25D0*(ONE-ETA)*(ONE+XI)
      P(1,5)=0.25D0*(ONE+ETA)*(ONE+XI)
      P(1,7)=0.25D0*(1.0D0*ONE+ETA)*(ONE-XI)
      P(2,2)=0.25D0*(1.0D0*ONE-ETA)*(1.0D0*ONE-XI)
      P(2,4)=0.25D0*(ONE-ETA)*(ONE+XI)
      P(2,6)=0.25D0*(ONE+ETA)*(ONE+XI)
      P(2,8)=0.25D0*(1.0D0*ONE+ETA)*(ONE-XI)      
C*****BJ:Jacobi Matrix
           DO I=1,2
             DO K=1,2
                BJ(I,K)=0.0D0
             ENDDO
           ENDDO
           
          CALL F_YAKEBI(BJ,BN,COORDS,xi,eta,NNODE-4)      !计算雅可比矩阵
          
C*****BR Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BR(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BH,BJ)
           
           BR(1,1)=BJ(2,2)/BH
           BR(1,2)=-1*BJ(1,2)/BH
           BR(2,1)=-1*BJ(2,1)/BH
           BR(2,2)=BJ(1,1)/BH
           
C           
C*****BL: Transition Matrix
           DO I=1,2
             DO K=1,4
               BL(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BL,BR,2,2,BN,2,4)  !BRxBN
C*****Make B the Matrix         
           DO I=1,3
             DO K=1,8
               B(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,4
             B(1,K*2-1)=BL(1,K)
             B(2,K*2)=BL(2,K)
             B(3,K*2-1)=BL(2,K)
             B(3,K*2)=BL(1,K)
           ENDDO  
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE SHAPE FUNCTION AND B Matrix OF VITURE ELEMENT   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
C*****BNV:Shape Functions     
           xiv=xi
           etav=eta
           DO I=1,2
           DO J=1,4
           BNV(I,J)=ZERO
           ENDDO
           ENDDO
C*****NV为形函数矩阵
            DO I=1,2
             DO J=1,8
               PV(I,J)=ZERO
             ENDDO
           ENDDO
      PV(1,1)=0.25D0*(1.0D0*ONE-ETAV)*(1.0D0*ONE-XIV)
      PV(1,3)=0.25D0*(ONE-ETAV)*(ONE+XIV)
      PV(1,5)=0.25D0*(ONE+ETAV)*(ONE+XIV)
      PV(1,7)=0.25D0*(1.0D0*ONE+ETAV)*(ONE-XIV)
      PV(2,2)=0.25D0*(1.0D0*ONE-ETAV)*(1.0D0*ONE-XIV)
      PV(2,4)=0.25D0*(ONE-ETAV)*(ONE+XIV)
      PV(2,6)=0.25D0*(ONE+ETAV)*(ONE+XIV)
      PV(2,8)=0.25D0*(1.0D0*ONE+ETAV)*(ONE-XIV)     
C*****BJV:Jacobi Matrix
           DO I=1,2
             DO K=1,2
                BJV(I,K)=0.0D0
             ENDDO
           ENDDO
           
          CALL F_YAKEBIV(BJV,BNV,COORDS,xiv,etav,4)      !计算雅可比矩阵
          
C*****BRV Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BRV(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BHV,BJV)
           
           BRV(1,1)=BJV(2,2)/BHV
           BRV(1,2)=-1*BJV(1,2)/BHV
           BRV(2,1)=-1*BJV(2,1)/BHV
           BRV(2,2)=BJV(1,1)/BHV
           
C           
C*****BLV: Transition Matrix
           DO I=1,2
             DO K=1,4
               BLV(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BLV,BRV,2,2,BNV,2,4)  !BRxBN
C*****Make BV the Matrix         
           DO I=1,3
             DO K=1,8
               BV(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,4
             BV(1,K*2-1)=BLV(1,K)
             BV(2,K*2)=BLV(2,K)
             BV(3,K*2-1)=BLV(2,K)
             BV(3,K*2)=BLV(1,K)
           END DO   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE COUPLE MATRIX CVJ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
      CALL F_TRANS(PVTRANSA,PV,2,8)
      CALL F_TRANS(BVTRANSA,BV,3,8)    
      CALL F_MATMUL(PVNJ,PVTRANSA,8,2,P,2,8)   
      CALL F_MATMUL(BVBJ,BVTRANSA,8,3,B,3,8)    
       DO I=1,8
            DO J=1,8
            CVJ(I,J)=(PVNJ(I,J)+BVBJ(I,J)*LLL*LLL)*BH
            END DO
       END DO     
      CALL F_TRANS(CVJTRANSA,CVJ,8,8)
           
             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE STRESSES AND STRAINS                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****INITIALIZE CDF_STRAIN,CDF_DSTRAIN, CDF_STRESS
           DO I=1,4
              CDF_STRAIN(I)=ZERO
              CDF_DSTRAIN(I)=ZERO
              CDF_STRESS(I)=ZERO
           ENDDO           
           DO K1=1,3
             DO K2=1,3
              D(K1,K2)=ZERO
             END DO
           END DO  	           
C
C*****由位移和位移增量求得应变和应变增量，同时把3X1向量，拓展到4X1
C     
           DO I=1,3
              DO J=1,8
                CDF_STRAIN(I)=CDF_STRAIN(I)+B(I,J)*U(J)
                CDF_DSTRAIN(I)=CDF_DSTRAIN(I)+B(I,J)*DU(J,1)
              ENDDO
           ENDDO 
C           
           CDF_STRAIN(4)=CDF_STRAIN(3)  
           CDF_STRAIN(3)=ZERO
           CDF_DSTRAIN(4)=CDF_DSTRAIN(3)  
           CDF_DSTRAIN(3)=ZERO
C
C*****从上一载荷步载入塑性应变、应力以及等效塑性应变的值
          
         DO K1=1,11
               UMAT_SVARS(K1)=SVARS(11*Kintk+K1-11)
          END DO                          	
C*****调用子函数F_UMAT材料矩阵以及输出应变、应力      
           CALL F_UMAT_MISES(CDF_STRAIN,CDF_DSTRAIN,NPROPS,PROPS,
     1            CDF_D,UMAT_SVARS)        
C
C***** MAKE S-MATRIX  
C  
C*****提取平面应变D材料矩阵，构建平面应变S阵                     
          DO I=1,3
                DO K=1,8
                   S(I,K)=0.0D0
              ENDDO
           ENDDO
           
           DO K1=1,2
             DO K2=1,2
              D(K1,K2)=CDF_D(K1,K2)
             END DO
              D(3,K1)=CDF_D(4,K1) 	
           END DO 
              D(3,3)=CDF_D(4,4)      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         MAKE S-MATRIX                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
           DO I=1,3
                DO K=1,8
                   DO J=1,3
                      S(I,K)=S(I,K)+D(I,J)*B(J,K)
                   ENDDO
              ENDDO
           ENDDO     
           
C*****提取应力构建平面应变应力向量（3X1）             
           DO I=1,2
             CDF_TSTRESS(I)=UMAT_SVARS(I+4)
           ENDDO 
             CDF_TSTRESS(3)=UMAT_SVARS(8)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE body force                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
         DO K1=1,8
		    BODYFORCE(K1,1)=0
         ENDDO
        	    BODYFORCE(2,1)=-P(1,1)*0*10*BH	
		    BODYFORCE(4,1)=-P(1,3)*0*10*BH
		    BODYFORCE(6,1)=-P(1,5)*0*10*BH
		    BODYFORCE(8,1)=-P(1,7)*0*10*BH	   
              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                     AFFIRM WEIGHT FUNCTION                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      XXJ=COORDS(1,1)*(1-xi)*(1-eta)/4+COORDS(1,2)*(1+xi)*(1-eta)/4+ 
     1 COORDS(1,3)*(1+xi)*(1+eta)/4+COORDS(1,4)*(1-xi)*(1+eta)/4
      YYJ=COORDS(2,1)*(1-xi)*(1-eta)/4+COORDS(2,2)*(1+xi)*(1-eta)/4+ 
     1 COORDS(2,3)*(1+xi)*(1+eta)/4+COORDS(2,4)*(1-xi)*(1+eta)/4
      
      weight=0.05      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                     ASSEMBLE RHS AND LHS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ASSEMBLE RHS AND LHS
     
           CALL F_TRANS(BTRANSA,B,3,8)
           CALL F_MATMUL(S1,BTRANSA,8,3,S,3,8) 
     
           DO K1=1,8
           RHS(K1,1)=RHS(K1,1)-(1-WEIGHT)*BODYFORCE(K1,1)
           END DO
  
            DO K1=1,8
              DO K2=1,8
                    AMATRX(K1,K2)=AMATRX(K1,K2)-(1-WEIGHT)*S1(K1,K2)*BH
                   ENDDO
               ENDDO  
        
              DO K1=9,16
              DO K2=1,8
                     AMATRX(K1,K2)=AMATRX(K1,K2)+CVJ(K1-8,K2)
                   ENDDO
                 ENDDO  
                 
              DO K1=1,8
              DO K2=9,16
                     AMATRX(K1,K2)=AMATRX(K1,K2)+CVJTRANSA(K1,K2-8)
                   ENDDO
                 ENDDO  
          
            DO K1=1,8            
            DO K2=1,3            
                  RHS(K1,1)=RHS(K1,1)+(1-WEIGHT)*
     1                  BTRANSA(K1,K2)*CDF_TSTRESS(K2)*BH 
            ENDDO   
            ENDDO 
     
            DO K1=1,8
            DO K2=1,8
                 RHS(K1,1)=RHS(K1,1)-CVJTRANSA(K1,K2)*U(K2+8)
            ENDDO  
             ENDDO
    
              DO K1=9,16
              DO K2=1,8
                    RHS(K1,1)=RHS(K1,1)-CVJ(K1-8,K2)*U(K2)
                   ENDDO
               ENDDO                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                  STORE STRAINS IN STATE VARIABLE ARRAY                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C***** 
C
           DO K1=1,11
               SVARS(11*Kintk+K1-11)=UMAT_SVARS(K1)
           END DO     
           
           DO I=1,11
              temporary_array(I)=UMAT_SVARS(I)*WEIGHT
           ENDDO
          
           DO I=1,11
              user_variable(jelem,I,Kintk)=temporary_array(I)
           END DO          
          
           IF(Kintk.EQ.1)THEN
           user_variable(jelem,12,Kintk)=COORDS(1,1)
           user_variable(jelem,13,Kintk)=COORDS(2,1)
           user_variable(jelem,14,Kintk)=(COORDS(1,2)+COORDS(1,1))/2
           user_variable(jelem,15,Kintk)=COORDS(2,2)  
           user_variable(jelem,16,Kintk)=(COORDS(1,1)+COORDS(1,2)
     1           +COORDS(1,3)+COORDS(1,4))/4
           user_variable(jelem,17,Kintk)=(COORDS(2,1)+COORDS(2,2)
     1           +COORDS(2,3)+COORDS(2,4))/4
           user_variable(jelem,18,Kintk)=(COORDS(1,4)+COORDS(1,1))/2
           user_variable(jelem,19,Kintk)=(COORDS(2,4)+COORDS(2,1))/2
           ELSE IF(Kintk.EQ.2)THEN
           user_variable(jelem,12,Kintk)=(COORDS(1,2)+COORDS(1,1))/2
           user_variable(jelem,13,Kintk)=COORDS(2,2) 
           user_variable(jelem,14,Kintk)=COORDS(1,2)  
           user_variable(jelem,15,Kintk)=COORDS(2,2)  
           user_variable(jelem,16,Kintk)=(COORDS(1,2)+COORDS(1,3))/2
           user_variable(jelem,17,Kintk)=(COORDS(2,2)+COORDS(2,3))/2
           user_variable(jelem,18,Kintk)=(COORDS(1,1)+COORDS(1,2)
     1           +COORDS(1,3)+COORDS(1,4))/4          
           user_variable(jelem,19,Kintk)=(COORDS(2,1)+COORDS(2,2)
     1           +COORDS(2,3)+COORDS(2,4))/4
                ELSE IF(Kintk.EQ.3)THEN
             user_variable(jelem,12,Kintk+1)=(COORDS(1,1)+COORDS(1,2)
     1           +COORDS(1,3)+COORDS(1,4))/4          
           user_variable(jelem,13,Kintk+1)=(COORDS(2,1)+COORDS(2,2)
     1           +COORDS(2,3)+COORDS(2,4))/4
           user_variable(jelem,14,Kintk+1)=(COORDS(1,2)+COORDS(1,3))/2
           user_variable(jelem,15,Kintk+1)=(COORDS(2,2)+COORDS(2,3))/2
           user_variable(jelem,16,Kintk+1)=COORDS(1,3)
           user_variable(jelem,17,Kintk+1)=COORDS(2,3)
           user_variable(jelem,18,Kintk+1)=(COORDS(1,4)+COORDS(1,3))/2
           user_variable(jelem,19,Kintk+1)=(COORDS(2,4)+COORDS(2,3))/2
            ELSE IF(Kintk.EQ.4)THEN
           user_variable(jelem,12,Kintk-1)=(COORDS(1,4)+COORDS(1,1))/2
           user_variable(jelem,13,Kintk-1)=(COORDS(2,4)+COORDS(2,1))/2
           user_variable(jelem,14,Kintk-1)=(COORDS(1,1)+COORDS(1,2)
     1           +COORDS(1,3)+COORDS(1,4))/4           
           user_variable(jelem,15,Kintk-1)=(COORDS(2,1)+COORDS(2,2)
     1           +COORDS(2,3)+COORDS(2,4))/4
           user_variable(jelem,16,Kintk-1)=(COORDS(1,4)+COORDS(1,3))/2
           user_variable(jelem,17,Kintk-1)=(COORDS(2,4)+COORDS(2,3))/2
           user_variable(jelem,18,Kintk-1)=COORDS(1,4)
           user_variable(jelem,19,Kintk-1)=COORDS(2,4)
           ENDIF          
      ENDDO
C*****guass point end  
      ENDIF
C      
C------------------------------------
C    精细单元与虚拟单元耦合
C----------------------------------      
       
      if (JTYPE.EQ.3) THEN    
      
      Gauss(1,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(1,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(2,1)=SQRT(ONE/THREE)
      Gauss(2,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(3,1)=SQRT(ONE/THREE)
      Gauss(3,2)=SQRT(ONE/THREE)
      Gauss(4,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(4,2)=SQRT(ONE/THREE)*1.0D0
C      
C------------------------------------
C    INITAIL  RHS、AMATRX
C    NDOFEL：节点的自由度数 
C    NRHS：  Number of load vectors.
C    AMATRX: 刚度矩阵
C    Gauss： 高斯积分
C----------------------------------
      DO K1 = 1, NDOFEL
           DO KRHS = 1, NRHS
             RHS(K1,KRHS) = ZERO
           ENDDO
           DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
           ENDDO
      ENDDO
     
       DO Kintk=1,4
           xi=Gauss(Kintk,1)
           eta=Gauss(Kintk,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Make B Matrix  in Natural Coordinate System                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****BN:Shape Functions     
           DO I=1,2
             DO J=1,4
               BN(I,J)=ZERO
             ENDDO
           ENDDO           
C*****N为形函数矩阵
            DO I=1,2
             DO J=1,8
               P(I,J)=ZERO
             ENDDO
           ENDDO
      P(1,1)=0.25D0*(1.0D0*ONE-ETA)*(1.0D0*ONE-XI)
      P(1,3)=0.25D0*(ONE-ETA)*(ONE+XI)
      P(1,5)=0.25D0*(ONE+ETA)*(ONE+XI)
      P(1,7)=0.25D0*(1.0D0*ONE+ETA)*(ONE-XI)
      P(2,2)=0.25D0*(1.0D0*ONE-ETA)*(1.0D0*ONE-XI)
      P(2,4)=0.25D0*(ONE-ETA)*(ONE+XI)
      P(2,6)=0.25D0*(ONE+ETA)*(ONE+XI)
      P(2,8)=0.25D0*(1.0D0*ONE+ETA)*(ONE-XI)      
C*****BJ:Jacobi Matrix
           DO I=1,2
             DO K=1,2
                BJ(I,K)=0.0D0
             ENDDO
           ENDDO
           
          CALL F_YAKEBI(BJ,BN,COORDS,xi,eta,NNODE-4)      !计算雅可比矩阵
          
C*****BR Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BR(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BH,BJ)
           
           BR(1,1)=BJ(2,2)/BH
           BR(1,2)=-1*BJ(1,2)/BH
           BR(2,1)=-1*BJ(2,1)/BH
           BR(2,2)=BJ(1,1)/BH
           
C           
C*****BL: Transition Matrix
           DO I=1,2
             DO K=1,4
               BL(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BL,BR,2,2,BN,2,4)  !BRxBN
C*****Make B the Matrix         
           DO I=1,3
             DO K=1,8
               B(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,4
             B(1,K*2-1)=BL(1,K)
             B(2,K*2)=BL(2,K)
             B(3,K*2-1)=BL(2,K)
             B(3,K*2)=BL(1,K)
           ENDDO  
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE SHAPE FUNCTION AND B Matrix OF VITURE ELEMENT   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
C*****BNV:Shape Functions     
       XXJ=COORDS(1,1)*(1-xi)*(1-eta)/4+COORDS(1,2)*(1+xi)*(1-eta)/4+ 
     1 COORDS(1,3)*(1+xi)*(1+eta)/4+COORDS(1,4)*(1-xi)*(1+eta)/4
        YYJ=COORDS(2,1)*(1-xi)*(1-eta)/4+COORDS(2,2)*(1+xi)*(1-eta)/4+ 
     1 COORDS(2,3)*(1+xi)*(1+eta)/4+COORDS(2,4)*(1-xi)*(1+eta)/4
        DO I=0,200
        DO J=0,200
        X1=-1+0.01*I
        Y1=-1+0.01*J
           F11=COORDS(1,5)*(1-X1)*(1-Y1)/4+COORDS(1,6)*(1+X1)*(1-Y1)/4+ 
     1 COORDS(1,7)*(1+X1)*(1+Y1)/4+COORDS(1,8)*(1-X1)*(1+Y1)/4-XXJ
           F12=COORDS(2,5)*(1-X1)*(1-Y1)/4+COORDS(2,6)*(1+X1)*(1-Y1)/4+ 
     1 COORDS(2,7)*(1+X1)*(1+Y1)/4+COORDS(2,8)*(1-X1)*(1+Y1)/4-YYJ      
      IF((ABS(F11).LE.0.001).AND.(ABS(F12).LE.0.001))THEN
      xiv=X1
      etav=Y1 
      GO TO 1000
      END IF
      END DO  
      END DO 
      
1000       DO I=1,2
           DO J=1,4
           BNV(I,J)=ZERO
           ENDDO
           ENDDO
C*****NV为形函数矩阵
           write(6,*)'xiv'
           write(6,*)(xiv)
           write(6,*)'etav'
           write(6,*)(etav)
           write(6,*)'jelem'
           write(6,*)(jelem)
           
            DO I=1,2
             DO J=1,8
               PV(I,J)=ZERO
             ENDDO
           ENDDO
      PV(1,1)=0.25D0*(1.0D0*ONE-ETAV)*(1.0D0*ONE-XIV)
      PV(1,3)=0.25D0*(ONE-ETAV)*(ONE+XIV)
      PV(1,5)=0.25D0*(ONE+ETAV)*(ONE+XIV)
      PV(1,7)=0.25D0*(1.0D0*ONE+ETAV)*(ONE-XIV)
      PV(2,2)=0.25D0*(1.0D0*ONE-ETAV)*(1.0D0*ONE-XIV)
      PV(2,4)=0.25D0*(ONE-ETAV)*(ONE+XIV)
      PV(2,6)=0.25D0*(ONE+ETAV)*(ONE+XIV)
      PV(2,8)=0.25D0*(1.0D0*ONE+ETAV)*(ONE-XIV)     
C*****BJV:Jacobi Matrix
           DO I=1,2
             DO K=1,2
                BJV(I,K)=0.0D0
             ENDDO
           ENDDO
           
          CALL F_YAKEBIV(BJV,BNV,COORDS,xiv,etav,4)      !计算雅可比矩阵
          
C*****BRV Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BRV(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BHV,BJV)
           
           BRV(1,1)=BJV(2,2)/BHV
           BRV(1,2)=-1*BJV(1,2)/BHV
           BRV(2,1)=-1*BJV(2,1)/BHV
           BRV(2,2)=BJV(1,1)/BHV
           
C           
C*****BLV: Transition Matrix
           DO I=1,2
             DO K=1,4
               BLV(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BLV,BRV,2,2,BNV,2,4)  !BRxBN
C*****Make BV the Matrix         
           DO I=1,3
             DO K=1,8
               BV(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,4
             BV(1,K*2-1)=BLV(1,K)
             BV(2,K*2)=BLV(2,K)
             BV(3,K*2-1)=BLV(2,K)
             BV(3,K*2)=BLV(1,K)
           END DO   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE COUPLE MATRIX CVJ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
      CALL F_TRANS(PVTRANSA,PV,2,8)
      CALL F_TRANS(BVTRANSA,BV,3,8)    
      CALL F_MATMUL(PVNJ,PVTRANSA,8,2,P,2,8)   
      CALL F_MATMUL(BVBJ,BVTRANSA,8,3,B,3,8)    
       DO I=1,8
            DO J=1,8
            CVJ(I,J)=(PVNJ(I,J)+BVBJ(I,J)*LLL*LLL)*BH
            END DO
       END DO     
      CALL F_TRANS(CVJTRANSA,CVJ,8,8)
           
             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE STRESSES AND STRAINS                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****INITIALIZE CDF_STRAIN,CDF_DSTRAIN, CDF_STRESS
           DO I=1,4
              CDF_STRAIN(I)=ZERO
              CDF_DSTRAIN(I)=ZERO
              CDF_STRESS(I)=ZERO
           ENDDO           
           DO K1=1,3
             DO K2=1,3
              D(K1,K2)=ZERO
             END DO
           END DO  	           
C
C*****由位移和位移增量求得应变和应变增量，同时把3X1向量，拓展到4X1
C     
           DO I=1,3
              DO J=1,8
                CDF_STRAIN(I)=CDF_STRAIN(I)+B(I,J)*U(J)
                CDF_DSTRAIN(I)=CDF_DSTRAIN(I)+B(I,J)*DU(J,1)
              ENDDO
           ENDDO 
C           
           CDF_STRAIN(4)=CDF_STRAIN(3)  
           CDF_STRAIN(3)=ZERO
           CDF_DSTRAIN(4)=CDF_DSTRAIN(3)  
           CDF_DSTRAIN(3)=ZERO
C
C*****从上一载荷步载入塑性应变、应力以及等效塑性应变的值
          
         DO K1=1,11
               UMAT_SVARS(K1)=SVARS(11*Kintk+K1-11)
          END DO                          	
C*****调用子函数F_UMAT材料矩阵以及输出应变、应力      
           CALL F_UMAT_MISES(CDF_STRAIN,CDF_DSTRAIN,NPROPS,PROPS,
     1            CDF_D,UMAT_SVARS)        
C
C***** MAKE S-MATRIX  
C  
C*****提取平面应变D材料矩阵，构建平面应变S阵                     
             DO I=1,3
                DO K=1,8
                   S(I,K)=0.0D0
              ENDDO
           ENDDO
           
           DO K1=1,2
             DO K2=1,2
              D(K1,K2)=CDF_D(K1,K2)
             END DO
              D(3,K1)=CDF_D(4,K1) 	
           END DO 
              D(3,3)=CDF_D(4,4)      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         MAKE S-MATRIX                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           
           DO I=1,3
                DO K=1,8
                   DO J=1,3
                      S(I,K)=S(I,K)+D(I,J)*B(J,K)
                   ENDDO
              ENDDO
           ENDDO  
           
C*****提取应力构建平面应变应力向量（3X1）             
           DO I=1,2
             CDF_TSTRESS(I)=UMAT_SVARS(I+4)
           ENDDO 
             CDF_TSTRESS(3)=UMAT_SVARS(8)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE body force                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
         DO K1=1,8
		    BODYFORCE(K1,1)=0
         ENDDO
        	    BODYFORCE(2,1)=-P(1,1)*0*10*BH	
		    BODYFORCE(4,1)=-P(1,3)*0*10*BH
		    BODYFORCE(6,1)=-P(1,5)*0*10*BH
		    BODYFORCE(8,1)=-P(1,7)*0*10*BH	   
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                     AFFIRM WEIGHT FUNCTION                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

            weight=0.95

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                     ASSEMBLE RHS AND LHS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ASSEMBLE RHS AND LHS
C
           CALL F_TRANS(BTRANSA,B,3,8)
           CALL F_MATMUL(S1,BTRANSA,8,3,S,3,8) 
          
           DO K1=1,8
           RHS(K1,1)=RHS(K1,1)-(1-WEIGHT)*BODYFORCE(K1,1)
           END DO
  
          DO K1=1,8
              DO K2=1,8
                    AMATRX(K1,K2)=AMATRX(K1,K2)-(1-WEIGHT)*S1(K1,K2)*BH
                   ENDDO
               ENDDO  
        
              DO K1=9,16
              DO K2=1,8
                     AMATRX(K1,K2)=AMATRX(K1,K2)-CVJ(K1-8,K2)
                   ENDDO
                 ENDDO  
                 
              DO K1=1,8
              DO K2=9,16
                     AMATRX(K1,K2)=AMATRX(K1,K2)-CVJTRANSA(K1,K2-8)
                   ENDDO
                 ENDDO  
      
       
            DO K1=1,8            
            DO K2=1,3            
          RHS(K1,1)=RHS(K1,1)+(1-WEIGHT)*BTRANSA(K1,K2)*CDF_TSTRESS(K2)
     1    *BH 
            ENDDO   
            ENDDO 
     
            DO K1=1,8
            DO K2=1,8
                 RHS(K1,1)=RHS(K1,1)+CVJTRANSA(K1,K2)*U(K2+8)
            ENDDO  
            ENDDO  
    
              DO K1=9,16
              DO K2=1,8
                    RHS(K1,1)=RHS(K1,1)+CVJ(K1-8,K2)*U(K2)
                   ENDDO
               ENDDO                   
   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                  STORE STRAINS IN STATE VARIABLE ARRAY                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C***** 
C
           DO K1=1,11
               SVARS(11*Kintk+K1-11)=UMAT_SVARS(K1)
           END DO     
           
           DO I=1,11
              temporary_array(I)=UMAT_SVARS(I)*WEIGHT
           ENDDO
          
           DO I=1,11
              user_variable(jelem,I,Kintk)=temporary_array(I)
           END DO          
          
           if(kintk.eq.1)then
           user_variable(jelem,12,Kintk)=XXJ
           user_variable(jelem,13,Kintk)=YYJ  
           endif
           
           if(kintk.eq.2)then
           user_variable(jelem,12,Kintk)=XXJ
           user_variable(jelem,13,Kintk)=YYJ  
           endif
           
           if(kintk.eq.3)then
           user_variable(jelem,12,Kintk+1)=XXJ
           user_variable(jelem,13,Kintk+1)=YYJ  
           endif
           
           if(kintk.eq.4)then
           user_variable(jelem,12,Kintk-1)=XXJ
           user_variable(jelem,13,Kintk-1)=YYJ  
           endif
           
      ENDDO
C*****guass point end  
      ENDIF
             
      if (JTYPE.EQ.4) THEN    
      
      Gauss(1,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(1,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(2,1)=SQRT(ONE/THREE)
      Gauss(2,2)=SQRT(ONE/THREE)*-1.0D0
      Gauss(3,1)=SQRT(ONE/THREE)
      Gauss(3,2)=SQRT(ONE/THREE)
      Gauss(4,1)=SQRT(ONE/THREE)*-1.0D0
      Gauss(4,2)=SQRT(ONE/THREE)*1.0D0
C      
C---------------------------------------
C  INITIALIZE  RHS、AMATRX
C  NDOFEL：节点的自由度数 
C  NRHS：  Number of load vectors.
C  AMATRX: 刚度矩阵
C  Gauss： 高斯积分
C---------------------------------------
      DO K1 = 1, NDOFEL
           DO KRHS = 1, NRHS
             RHS(K1,KRHS) = ZERO
           ENDDO
           DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
           ENDDO
      ENDDO
      
      DO Kintk=1,NNODE
           xi=Gauss(Kintk,1)
           eta=Gauss(Kintk,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Make B Matrix  in Global Coordinate System                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*****BN:Shape Functions     
           DO I=1,2
             DO J=1,NNODE
               BN(I,J)=ZERO
             ENDDO
           ENDDO

C*****BJ:Jacobi Matrix
          DO I=1,2
             DO J=1,8
               P(I,J)=ZERO
             ENDDO
           ENDDO
      P(1,1)=0.25*(1.0*ONE-ETA)*(1.0*ONE-XI)
      P(1,3)=0.25*(ONE-ETA)*(ONE+XI)
      P(1,5)=0.25*(ONE+ETA)*(ONE+XI)
      P(1,7)=0.25*(1.0*ONE+ETA)*(ONE-XI)
      P(2,2)=0.25*(1.0*ONE-ETA)*(1.0*ONE-XI)
      P(2,4)=0.25*(ONE-ETA)*(ONE+XI)
      P(2,6)=0.25*(ONE+ETA)*(ONE+XI)
      P(2,8)=0.25*(1.0*ONE+ETA)*(ONE-XI)  
         DO I=1,2
             DO K=1,2
                BJ(I,K)=0.0D0
             ENDDO
           ENDDO
           
           CALL F_YAKEBI(BJ,BN,COORDS,xi,eta,4)   !计算雅可比矩阵
          
C*****BR Jacobi Inverse Matrix   
           DO I=1,2
              DO K=1,2
                   BR(I,K)=0.0D0
              ENDDO
           ENDDO
C           
           CALL F_DETA(BH,BJ)
           
           BR(1,1)=BJ(2,2)/BH
           BR(1,2)=-1*BJ(1,2)/BH
           BR(2,1)=-1*BJ(2,1)/BH
           BR(2,2)=BJ(1,1)/BH
           
C           
C*****BL: Transition Matrix
           DO I=1,2
             DO K=1,NNODE
               BL(I,K)=0.0D0
             ENDDO
           ENDDO 
           CALL F_MATMUL(BL,BR,2,2,BN,2,4)  !BRxBN
C*****Make B the Matrix         
           DO I=1,3
             DO K=1,8
               B(I,K)=0.0D0
             ENDDO
           ENDDO
               
           DO K=1,NNODE
             B(1,K*2-1)=BL(1,K)
             B(2,K*2)=BL(2,K)
             B(3,K*2-1)=BL(2,K)
             B(3,K*2)=BL(1,K)
           ENDDO           
   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE STRESSES AND STRAINS                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*****INITIALIZE CDF_STRAIN,CDF_DSTRAIN, CDF_STRESS
           DO I=1,4
              CDF_STRAIN(I)=ZERO
              CDF_DSTRAIN(I)=ZERO
              CDF_STRESS(I)=ZERO
           ENDDO

           
           DO K1=1,3
             DO K2=1,3
              D(K1,K2)=ZERO
             END DO
           END DO  	           
C
C*****由位移和位移增量求得应变和应变增量，同时把3X1向量，拓展到4X1
C     
           DO I=1,3
              DO J=1,8
                CDF_STRAIN(I)=CDF_STRAIN(I)+B(I,J)*U(J)
                CDF_DSTRAIN(I)=CDF_DSTRAIN(I)+B(I,J)*DU(J,1)
              ENDDO
           ENDDO 
C           
           CDF_STRAIN(4)=CDF_STRAIN(3)  
           CDF_STRAIN(3)=ZERO
           CDF_DSTRAIN(4)=CDF_DSTRAIN(3)  
           CDF_DSTRAIN(3)=ZERO
C
C*****从上一载荷步载入塑性应变、应力以及等效塑性应变的值
C  
           DO K1=1,11
               UMAT_SVARS(K1)=SVARS(11*Kintk+K1-11)
           END DO  
           	
C*****调用子函数F_UMAT材料矩阵以及输出塑性应变、流动应力  
C*****和等效塑性应变          
           CALL F_UMAT_MISES(CDF_STRAIN,CDF_DSTRAIN,NPROPS,PROPS,
     +            CDF_D,UMAT_SVARS)        
C
C***** MAKE S-MATRIX  
C  
C*****提取平面应变D材料矩阵，构建平面应变S阵                     
            DO I=1,3
                DO K=1,8
                   S(I,K)=0.0D0
              ENDDO
           ENDDO
           
           
           DO K1=1,2
             DO K2=1,2
              D(K1,K2)=CDF_D(K1,K2)
             END DO
              D(3,K1)=CDF_D(4,K1)
              D(3,K1)=CDF_D(4,K1)             	
           END DO 
	      D(3,3)=CDF_D(4,4)      
          DO I=1,3
                DO K=1,8
                   DO J=1,3
                      S(I,K)=S(I,K)+D(I,J)*B(J,K)
                   ENDDO
              ENDDO
           ENDDO       
C*****提取应力构建平面应变应力向量（3X1）             
           DO I=1,2
             CDF_TSTRESS(I)=UMAT_SVARS(I+4)
           ENDDO 
             CDF_TSTRESS(3)=UMAT_SVARS(8)                     
CC         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               CALCULATE body force                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		   
         DO K1=1,8
		    BODYFORCE(K1,1)=0
         ENDDO
        	    BODYFORCE(2,1)=-P(1,1)*2.5*10*BH	
		    BODYFORCE(4,1)=-P(1,3)*2.5*10*BH
		    BODYFORCE(6,1)=-P(1,5)*2.5*10*BH
		    BODYFORCE(8,1)=-P(1,7)*2.5*10*BH	  
C*****ASSEMBLE RHS AND LHS   
C                       
C*****组装RHS和AMATRX
          
          CALL F_TRANS(BTRANSA,B,3,8)
            
          DO K1=1,8
          RHS(K1,1)=RHS(K1,1)+BODYFORCE(K1,1)
          END DO

           DO K1=1,8
			    DO KK=1,3
                 RHS(K1,1)=RHS(K1,1)-BTRANSA(K1,KK)*CDF_TSTRESS(KK)*BH
                 ENDDO
          ENDDO
          
          DO K1=1,8
              DO K2=1,8
                 DO KK=1,3
              AMATRX(K1,K2)=AMATRX(K1,K2)+BTRANSA(K1,KK)*S(KK,K2)*BH
                 ENDDO
              ENDDO
           ENDDO
C                      
C***** STORE STRAINS IN STATE VARIABLE ARRAY
C
           DO K1=1,11
               SVARS(11*Kintk+K1-11)=UMAT_SVARS(K1)
           END DO     
           
           DO I=1,11
              temporary_array(I)=UMAT_SVARS(I)  
           ENDDO    
             
             
          if(kintk.eq.1)then
           do i=1,11
           user_variable(jelem,I,Kintk)=temporary_array(I)
           enddo
           endif
           
           if(kintk.eq.2)then
            do i=1,11
           user_variable(jelem,I,Kintk)=temporary_array(I)
           enddo
           endif
    
           if(kintk.eq.3)then
             do i=1,11
           user_variable(jelem,I,Kintk+1)=temporary_array(I)
           enddo
           endif
           
           if(kintk.eq.4)then
             do i=1,11
           user_variable(jelem,I,Kintk-1)=temporary_array(I)
           enddo
           endif                 
   
      ENDDO
C*****guass point end  
      ENDIF     
      
      RETURN
      END    

C##################################################################################################
C                                                                                                 #
C                                    SUBROUTINE LIST                                              #
C                                                                                                 #
C                                                                                                 #
C##################################################################################################
C 
C-----------------------------------------------------------------------------C      
C**************************邓肯张程序***************************
C-----------------------------------------------------------------------------C           
      SUBROUTINE F_UMAT(UMAT_STRAN,UMAT_DSTRAN,NPROPS,UMAT_PROPS,
     +            UMAT_D,UMAT_SVARS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER (ZERO= 0.D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     + NEWTON=10,TOLER=1.D-6)  
      DIMENSION UMAT_STRESS(4),UMAT_STRAN(4),UMAT_DSTRAN(4),
     +  UMAT_D(4,4),UMAT_EELAS(4),UMAT_PROPS(*),UMAT_EPLAS(4),
     +  UMAT_FLOW(4),UMAT_SVARS(11),UMAT_DSTRESS(4),ps(3)             
C -----------------------------------------------------------
C     UMAT FOR ISOTROPIC ELASTICITY AND ISOTROPIC PLASTICITY
C     PROPS(1) - E
C     PROPS(2) - NU
C     PROPS(3) - SYIELD
C     CALLS AHARD FOR CURVE OF SYIELD VS. PEEQ
C -----------------------------------------------------------
C
      ek=UMAT_PROPS(1)
	en=UMAT_PROPS(2)
	rf=UMAT_PROPS(3)
	c=UMAT_PROPS(4)
	phi=UMAT_PROPS(5)/180.0*3.1415926
	ug=UMAT_PROPS(6)
	ud=UMAT_PROPS(7)
	uf=UMAT_PROPS(8)
	ekur=UMAT_PROPS(9)
	pa=UMAT_PROPS(10)
	dphi=UMAT_PROPS(11)/180.0*3.1415926
      S1S3O=UMAT_SVARS(9)     !largest value of loading function
	S3O=UMAT_SVARS(10)      !consolidation stress
	SSS=UMAT_SVARS(11)       !stress level
   
      if(S30.lt.1.0)then
      S30=100.0
      endif
      
      DO K1=1,4
      UMAT_STRESS(K1)=UMAT_SVARS(K1+4)
      ENDDO
                
      call getps(UMAT_STRESS,ps,4)      
      phi=phi-dphi*log10(S3O/pa)
      
      call getemod(ps,ek,en,rf,c,phi,enu,pa,ekur,emod,s,S3O,ug,ud,uf
     1	,SSS,S1S3O)	     
  
      ebulk3=emod/(one-two*enu)
      eg2=emod/(one+enu)
      eg=eg2/two
      eg3=three*eg
      elam=(ebulk3-eg2)/three
  
    
 
      call getddsdde(UMAT_D,4,3,elam,eg2,eg)
      UMAT_DSTRESS=0.0
	call getstress(UMAT_D,UMAT_DSTRESS,UMAT_DSTRAN,4)
      do 701 i1=1,4
	UMAT_STRESS(i1)=UMAT_STRESS(i1)+UMAT_DSTRESS(i1)
701	continue
      
      
      call getps(UMAT_STRESS,ps,4)
      call getemod(ps,ek,en,rf,c,phi,enu,pa,ekur,emod,s,S3O,ug,ud,uf,
     1	SSS,S1S3O)
      ebulk3=emod/(one-two*enu)
      eg2=emod/(one+enu)
      eg=eg2/two
      eg3=three*eg
      elam=(ebulk3-eg2)/three
	call getddsdde(UMAT_D,4,3,elam,eg2,eg)
      if(ps(3).gt.S3O)S3O=ps(3)
      IF((PS(1)-PS(3)).GT.S1S3O)S1S3O=PS(1)-PS(3)
	IF(S.GT.SSS)SSS=S
      if(sss.LT.0.05) sss=0.05
      if((sss-0.95).GT.0) sss=0.95
      UMAT_SVARS(9)=S1S3O     !largest value of loading function
	UMAT_SVARS(10)=S3O     !consolidation stress
	UMAT_SVARS(11)=SSS     !stress level
      
      DO K1=1,4
        UMAT_SVARS(K1)=UMAT_STRAN(K1)
        UMAT_SVARS(K1+4)=UMAT_STRESS(K1)
      ENDDO 
      
      RETURN
      END
 
      subroutine getps(stress,ps,ntens) 
	include 'aba_param.inc'
	dimension ps(3),stress(ntens)    
      P=(STRESS(2)-STRESS(1))/2.0				 !求主应力
      Q=STRESS(1)+P
      R=SQRT(STRESS(4)*STRESS(4)+P*P)
      sigma01=Q+R
      sigma03=Q-R   
      if (sigma03.GT.sigma01) then 
      temp=sigma03
      sigma03=sigma01
      sigma01=temp
      endif      
      ps(1)=sigma03
      ps(2)=0
      ps(3)=sigma01      
	do 330 k1=1,3
	ps(k1)=-ps(k1)
330	continue	
      return
	end 

	subroutine getemod(ps,ek,en,rf,c,phi,enu,pa,ekur,emod,s,S3O
     1	,ug,ud,uf,SSS,S1S3O)
	include 'aba_param.inc'
	dimension ps(3)
	s=(1-sin(phi))*(ps(1)-ps(3))
	if(ps(3).lt.0) then
	psfei=0.1
	else
	psfei=ps(3)
	end if
	s=s/(2*c*cos(phi)+2*psfei*sin(phi))

      if(s.ge.0.99) then
	s=0.99
	end if
	
      aa=ud*(ps(1)-ps(3))
	aa=aa/(ek*pa*((S3O/pa)**en))
	aa=aa/(1-rf*s)
	enu=ug-uf*log10(S3O/pa)
	enu=enu/(1-aa)/(1-aa)
      if(enu.gt.0.49)enu=0.49
	if(enu.lt.0.01)enu=0.01
	
      EMOD=EK*PA*((S3O/PA)**EN)*((1-RF*S)**2)
  
	IF((S.LT.SSS).AND.((PS(1)-PS(3)).LT.S1S3O))THEN
	EMOD=EKUR*PA*((S3O/PA)**EN)
	END IF
     
      return
      end 

	subroutine getddsdde(ddsdde,ntens,ndi,elam,eg2,eg)
	include 'aba_param.inc'
	dimension ddsdde(ntens,ntens)
	do 20 k1=1,ntens
        do 10 k2=1,ntens
           ddsdde(k2,k1)=0.0
 10     continue
 20   continue
      do 40 k1=1,ndi
        do 30 k2=1,ndi
           ddsdde(k2,k1)=elam
 30     continue
        ddsdde(k1,k1)=eg2+elam
 40   continue
      do 50 k1=ndi+1,ntens
      ddsdde(k1,k1)=eg
 50   continue	
      return
	end   
      
      subroutine getstress(ddsdde,stress,dstran,ntens)
	include 'aba_param.inc'
	dimension ddsdde(ntens,ntens),stress(ntens),dstran(ntens)
	do 70 k1=1,ntens
        do 60 k2=1,ntens
           stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1)
 60     continue
 70   continue
	return
	end
      
C-----------------------------------------------------------------------------C      
C    用户自定义材料子程序，各向同性硬化，弹塑性，J2流动理论，适用于平面应变
C    IN：应变UMAT_STRAN(4),UMAT_DSTRAN(4)
C       
C    
C    OUT：材料矩阵
C
C    更新等效塑性应变
C-----------------------------------------------------------------------------C           
      SUBROUTINE F_UMAT_MISES(UMAT_STRAN,UMAT_DSTRAN,NPROPS,UMAT_PROPS,
     +            UMAT_D,UMAT_SVARS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER (ZERO= 0.D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     + NEWTON=10,TOLER=1.D-6)  
      DIMENSION UMAT_STRESS(4),UMAT_STRAN(4),UMAT_DSTRAN(4),
     +  UMAT_D(4,4),UMAT_EELAS(4),UMAT_PROPS(*),UMAT_EPLAS(4),
     +  UMAT_FLOW(4),UMAT_SVARS(9)               
C -----------------------------------------------------------
C     UMAT FOR ISOTROPIC ELASTICITY AND ISOTROPIC PLASTICITY
C     PROPS(1) - E
C     PROPS(2) - NU
C     PROPS(3) - SYIELD
C     CALLS AHARD FOR CURVE OF SYIELD VS. PEEQ
C -----------------------------------------------------------
C
      EMOD=UMAT_PROPS(1)
      ENU=UMAT_PROPS(2)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C     ELASTIC STIFFNESS
C
      DO K1=1,4
        DO K2=1,4
           UMAT_D(K2,K1)=0.0
        ENDDO
      ENDDO
C
      DO K1=1,3
        DO K2=1,3
           UMAT_D(K2,K1)=ELAM
        ENDDO
           UMAT_D(K1,K1)=EG2+ELAM
      ENDDO
           UMAT_D(4,4)=EG
C
C    CALCULATE STRESS FROM ELASTIC STRAINS
C
      DO K1=1,4
        UMAT_STRESS(K1)=UMAT_SVARS(K1+4)
      ENDDO
      DO K1=1,4
        DO K2=1,4
         UMAT_STRESS(K1)=UMAT_STRESS(K1)+UMAT_D(K1,K2)*UMAT_DSTRAN(K2)
        ENDDO
      ENDDO  
C
C    RECOVER ELASTIC AND PLASTIC STRAINS
C
      DO K1=1,4
         UMAT_EPLAS(K1)=UMAT_SVARS(K1)
         UMAT_EELAS(K1)=UMAT_STRAN(K1)-UMAT_SVARS(K1)+UMAT_DSTRAN(K1)   
      ENDDO
         EQPLAS=UMAT_SVARS(1+2*4)
C
C    IF NO YIELD STRESS IS GIVEN, MATERIAL IS TAKEN TO BE ELASTIC
C
      IF(NPROPS.GT.2.AND.UMAT_PROPS(3).GT.0.0) THEN
C
C       MISES STRESS
C
        SMISES=(UMAT_STRESS(1)-UMAT_STRESS(2))*(UMAT_STRESS(1)-
     +          UMAT_STRESS(2))+(UMAT_STRESS(2)-UMAT_STRESS(3))*
     +         (UMAT_STRESS(2)-UMAT_STRESS(3))+(UMAT_STRESS(3)-
     +          UMAT_STRESS(1))*(UMAT_STRESS(3)-UMAT_STRESS(1))    
        
        SMISES=SMISES+SIX*UMAT_STRESS(4)*UMAT_STRESS(4)
 
        SMISES=SQRT(SMISES/TWO)
C
C       HARDENING CURVE, GET YIELD STRESS
C
        NVALUE=NPROPS/2-1
        SYIEL0=0
        CALL F_HARDSUB(SYIEL0,HARD,EQPLAS,UMAT_PROPS(3),NVALUE)
C
C       DETERMINE IF ACTIVELY YIELDING
C
        IF (SMISES.GT.(1.0+TOLER)*SYIEL0) THEN
C
C         FLOW DIRECTION
C
          SHYDRO=(UMAT_STRESS(1)+UMAT_STRESS(2)+UMAT_STRESS(3))/THREE
          ONESY=ONE/SMISES
          DO K1=1,3
             UMAT_FLOW(K1)=ONESY*(UMAT_STRESS(K1)-SHYDRO)
          ENDDO
          DO K1=4,4
             UMAT_FLOW(K1)=UMAT_STRESS(K1)*ONESY
          ENDDO
C
C       SOLVE FOR EQUIV STRESS, NEWTON ITERATION
C
          SYIELD=SYIEL0
          DEQPL=0.0
          DO KEWTON=1,NEWTON
             RHS=SMISES-EG3*DEQPL-SYIELD
             DEQPL=DEQPL+RHS/(EG3+HARD)
             CALL F_HARDSUB(SYIELD,HARD,EQPLAS+DEQPL,UMAT_PROPS(3),
     +                    NVALUE)
             IF(ABS(RHS).LT.TOLER*SYIEL0) GOTO 103
          ENDDO
          WRITE(6,204) NEWTON
 204      FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     +        'CONVERGE AFTER ',I3,' ITERATIONS')
 103      CONTINUE
          EFFHRD=EG3*HARD/(EG3+HARD)
C
C       CALC STRESS AND UPDATE STRAINS
C
          DO K1=1,3
           UMAT_STRESS(K1)=UMAT_FLOW(K1)*SYIELD+SHYDRO
           UMAT_EPLAS(K1)=UMAT_EPLAS(K1)+THREE*UMAT_FLOW(K1)*DEQPL/TWO
           UMAT_EELAS(K1)=UMAT_EELAS(K1)-THREE*UMAT_FLOW(K1)*DEQPL/TWO
          ENDDO
          DO K1=3+1,4
             UMAT_STRESS(K1)=UMAT_FLOW(K1)*SYIELD
             UMAT_EPLAS(K1)=UMAT_EPLAS(K1)+THREE*UMAT_FLOW(K1)*DEQPL
             UMAT_EELAS(K1)=UMAT_EELAS(K1)-THREE*UMAT_FLOW(K1)*DEQPL
          ENDDO
          EQPLAS=EQPLAS+DEQPL
          SPD=DEQPL*(SYIEL0+SYIELD)/TWO
C
C       JACOBIAN
C
          EFFG=EG*SYIELD/SMISES
          EFFG2=TWO*EFFG
          EFFG3=THREE*EFFG2/TWO
          EFFLAM=(EBULK3-EFFG2)/THREE
          DO K1=1,3
             DO K2=1,3
                UMAT_D(K2,K1)=EFFLAM
             ENDDO
             UMAT_D(K1,K1)=EFFG2+EFFLAM
          ENDDO 
          UMAT_D(4,4)=EFFG

          DO K1=1,4
             DO  K2=1,4
                UMAT_D(K2,K1)=UMAT_D(K2,K1)+UMAT_FLOW(K2)
     +                      *UMAT_FLOW(K1)*(EFFHRD-EFFG3)
              ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C STORE STRAINS IN STATE VARIABLE ARRAY
C
      DO K1=1,4
        UMAT_SVARS(K1)=UMAT_EPLAS(K1)
        UMAT_SVARS(K1+4)=UMAT_STRESS(K1)
      ENDDO 
        UMAT_SVARS(1+2*4)=EQPLAS
C  
      RETURN
      END
C-----------------------------------------------------------------------------C      
c  子程序，根据等效塑性应变，利用插值的方法得到对应的屈服应力
C-----------------------------------------------------------------------------C
      SUBROUTINE F_HARDSUB(SYIELD,HARD,EQPLAS,TABLE,NVALUE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER (ZERO= 0.D0)   
      DIMENSION TABLE(2,NVALUE)
C
C
C    SET YIELD STRESS TO LAST VALUE OF TABLE, HARDENING TO ZERO
C
      
      SYIELD=TABLE(1,NVALUE)
      HARD=ZERO
C
C    IF MORE THAN ONE ENTRY, SEARCH TABLE
C
      IF(NVALUE.GT.1) THEN
        DO 101 K1=1,NVALUE-1
           EQPL1=TABLE(2,K1+1)
           IF(EQPLAS.LT.EQPL1) THEN
             EQPL0=TABLE(2,K1)
             IF(EQPL1.LE.EQPL0) THEN
                WRITE(6,203)
203             FORMAT(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     +                 'ENTERED IN ASCENDING ORDER')
                CALL XIT
             ENDIF
              DEQPL=EQPL1-EQPL0
              SYIEL0=TABLE(1,K1)
              SYIEL1=TABLE(1,K1+1)
              DSYIEL=SYIEL1-SYIEL0
              HARD=DSYIEL/DEQPL
              SYIELD=SYIEL0+(EQPLAS-EQPL0)*HARD
              GOTO 102
           ENDIF
101     CONTINUE
102     CONTINUE
      ENDIF
      RETURN
      END      
      
C-------------------------------------------------------------------------------------------------C
C  计算雅可比行列式
C  USER SUBROUTINE YAKEBI
C   IN : XKSI YETA (高斯积分点在母单元的坐标) XY (节点坐标)(2,4)  MNODE :节点数
C   OUT: YAKEBII,BN
C--------------------------------------------------------------------------------------------------C
      SUBROUTINE F_YAKEBI(YAKEBI,BSF,XY,XKSI,YETA,MNODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (KCOOR=2,ONE = 1.D0)
      DIMENSION YAKEBI(2,2),XY(2,4)
      DIMENSION BSF(2,4)
C*****BSF 为N矩阵
      BSF(1,1)=0.25D0*(-1.0D0*ONE+YETA)
      BSF(1,2)=0.25D0*(ONE-YETA)
      BSF(1,3)=0.25D0*(ONE+YETA)
      BSF(1,4)=0.25D0*(-1.0D0*ONE-YETA)
      BSF(2,1)=0.25D0*(-1.0D0*ONE+XKSI)
      BSF(2,2)=0.25D0*(-1.0D0*ONE-XKSI)
      BSF(2,3)=0.25D0*(ONE+XKSI)
      BSF(2,4)=0.25D0*(ONE-XKSI)
C*****YAKEBIIan Matrix
      DO J=1,MNODE
        YAKEBI(1,1)=YAKEBI(1,1)+BSF(1,J)*XY(1,J)
        YAKEBI(1,2)=YAKEBI(1,2)+BSF(1,J)*XY(2,J)
        YAKEBI(2,1)=YAKEBI(2,1)+BSF(2,J)*XY(1,J)
        YAKEBI(2,2)=YAKEBI(2,2)+BSF(2,J)*XY(2,J)
      ENDDO
      RETURN
      END 
C       
C-------------------------------------------------------------------------------------------------C
C  计算虚拟单元雅可比行列式
C  USER SUBROUTINE YAKEBIV
C   IN : XKSI YETA (高斯积分点在母单元的坐标) XY (节点坐标)(2,4)  MNODE :节点数
C   OUT: YAKEBII,BN
C--------------------------------------------------------------------------------------------------C
      SUBROUTINE F_YAKEBIV(YAKEBI,BSF,XY,XKSI,YETA,MNODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (KCOOR=2,ONE = 1.D0)
      DIMENSION YAKEBI(2,2),XY(2,8)
      DIMENSION BSF(2,4)
    
C*****BSF 为N矩阵
      BSF(1,1)=0.25D0*(-1.0D0*ONE+YETA)
      BSF(1,2)=0.25D0*(ONE-YETA)
      BSF(1,3)=0.25D0*(ONE+YETA)
      BSF(1,4)=0.25D0*(-1.0D0*ONE-YETA)
      BSF(2,1)=0.25D0*(-1.0D0*ONE+XKSI)
      BSF(2,2)=0.25D0*(-1.0D0*ONE-XKSI)
      BSF(2,3)=0.25D0*(ONE+XKSI)
      BSF(2,4)=0.25D0*(ONE-XKSI)
C*****YAKEBIIan Matrix
      DO J=1,MNODE
        YAKEBI(1,1)=YAKEBI(1,1)+BSF(1,J)*XY(1,J+4)
        YAKEBI(1,2)=YAKEBI(1,2)+BSF(1,J)*XY(2,J+4)
        YAKEBI(2,1)=YAKEBI(2,1)+BSF(2,J)*XY(1,J+4)
        YAKEBI(2,2)=YAKEBI(2,2)+BSF(2,J)*XY(2,J+4)
      ENDDO
      RETURN
      END 
C       
C-------------------------------------------------------------------------------------------------C
C   矩阵相乘
C   USER SUBROUTINE F_MATMUL
C   IN : XKSI YETA (高斯积分点在母单元的坐标) XY (节点坐标)(2,4)  MNODE :节点数
C   OUT: YAKEBII,BN
C--------------------------------------------------------------------------------------------------C
C       矩阵乘法的子程序
C
      SUBROUTINE F_MATMUL(AMATRIX,BMATRIX,MBROW,MBCOL,CMATRIX,
     1     MCROW,MCCOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BMATRIX(MBROW,MBCOL),CMATRIX(MCROW,MCCOL)
      DIMENSION AMATRIX(MBROW,MCCOL)
C
      IF ( MBCOL .NE. MCROW ) THEN
        WRITE(*,*) 'Matrix size error!'
        STOP
      END IF
 
      DO I=1,MBROW
        DO J=1,MCCOL
          AMATRIX(I,J)=0
          DO K=1,MBCOL
            AMATRIX(I,J)=AMATRIX(I,J)+BMATRIX(I,K)*CMATRIX(K,J)
          END DO
        END DO
      END DO
 
      RETURN
      END
C-------------------------------------------------------------------------------------------------C
C   USER SUBROUTINE F_DETA:2X2 矩阵值的计算
C--------------------------------------------------------
C   Objective: 计算det[A]
C   IN: MATRIX	AMATRIX(2,2)
C   OUT: DETERMINATION OF MATERIXA -- ADET
C-------------------------------------------------------------------------------------------------C
      SUBROUTINE F_DETA(ADET,AMATRIX)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
      DIMENSION AMATRIX(2,2)
      TEMP1 = AMATRIX(1,1)*AMATRIX(2,2)
      TEMP2 = AMATRIX(1,2)*AMATRIX(2,1)
      ADET = TEMP1-TEMP2
      RETURN
      END   
      
C--------------------------------------------------------------------------------------------------C
C   USER SUBROUTINE TRANS: 矩阵转置
C   IN : A,NAROW,NACOL
C   out: TRANSA
C--------------------------------------------------------------------------------------------------C
      SUBROUTINE F_TRANS(TRANSA,AMATRIX,NAROW,NACOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
      DIMENSION AMATRIX(NAROW,NACOL),TRANSA(NACOL,NAROW)
      DO I = 1,NAROW
	 DO J = 1,NACOL
	    TRANSA(J,I) = AMATRIX(I,J)
	 ENDDO
      ENDDO
      RETURN
      END 	                 
          
C##############################################################################       
C     uel postprocess
C##############################################################################      
       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1                 NUVARM,NOEL,NPT,NLAYER,NSPT,KSTEP,KINC,
     2                 NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO, LACCFLG)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      DIMENSION UVAR(*),TIME(2),DIRECT(3,3),T(3,3),COORD(*),
     $     JMAC(*),JMATYP(*) 
C     USER DEFINED DIMENSION STATEMENTS
      CHARACTER*3 FLGRAY(15)
      DIMENSION ARRAY(15),JARRAY(15)
      PARAMETER (nelem=1000000,ninpt=4,nsvint=19)
      common /kuser/user_variable(nelem,nsvint,ninpt)
      
  
      IF(noel.GE.2289.AND.noel.LE.2344)THEN
      knoel=noel-2288
      ELSE IF(noel.GE.2345.AND.noel.LE.2360)THEN
      knoel=noel-472
      ELSE IF(noel.GE.2361.AND.noel.LE.3760)THEN
      knoel=noel-2288
      ELSE IF(noel.GE.3761.AND.noel.LE.4160)THEN
      knoel=noel-1872
      END IF
   
c      knoel=noel-2288

      DO I=1,19
      uvar(I)=user_variable(knoel,I,NPT)
      end do  
   
         IF(knoel.GE.1889.AND.knoel.LE.2288)THEN
               DO I=1873,1888
              DO K=1,4
           IF(user_variable(knoel,13,NPT).GE.user_variable(I,13,K).AND.
     1       user_variable(knoel,13,NPT).LE.user_variable(I,19,K).AND.
     1    (user_variable(knoel,12,NPT).GE.user_variable(I,12,K).OR.
     1    (user_variable(knoel,12,NPT).LT.user_variable(I,12,K).AND.
     1    user_variable(knoel,12,NPT).GE.user_variable(I,18,K)))
     1    .AND.(user_variable(knoel,12,NPT).LE.user_variable(I,14,K).OR.
     1     (user_variable(knoel,12,NPT).GT.user_variable(I,14,K).AND.
     1          user_variable(knoel,12,NPT).LE.user_variable(I,16,K)))
     1       )THEN
          DO J=5,9
         uvar(J)=user_variable(knoel,J,NPT)+user_variable(I,J,K)
          ENDDO
        ENDIF   
         ENDDO
            ENDDO
      ENDIF

      call getpps(uvar,19)    
       SMISES=(uvar(5)-uvar(6))*(uvar(5)-
     +          uvar(6))+(uvar(6)-uvar(7))*
     +         (uvar(6)-uvar(7))+(uvar(7)-
     +          uvar(5))*(uvar(7)-uvar(5))           
      SMISES=SMISES+6.0*uvar(8)*uvar(8)
      SMISES=SQRT(SMISES/2.0)    
      uvar(12)=SMISES  
      
      return
      end
      
      subroutine getpps(stress,ntens) 
	include 'aba_param.inc'
	dimension ps(3),stress(ntens)    
      P=(STRESS(6)-STRESS(5))/2.0				 !求主应力
      Q=STRESS(5)+P
      R=SQRT(STRESS(8)*STRESS(8)+P*P)
      sigma01=Q+R
      sigma03=Q-R   
      if (sigma03.GT.sigma01) then  
      temp=sigma03
      sigma03=sigma01
      sigma01=temp      
      endif         
      stress(10)=sigma01
      stress(11)=sigma03      
      return
	end 
C##############################################################################       
C     surface traction
C##############################################################################           
      
      SUBROUTINE UTRACLOAD(ALPHA,T_USER,KSTEP,KINC,TIME,NOEL,NPT,
     1 COORDS,DIRCOS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION T_USER(3), TIME(2), COORDS(3), DIRCOS(3,3)
      CHARACTER*80 SNAME

c      ALPHA=1/2/(1*1*1/12)*(1*1/4-coords(2)*coords(2))


      ALPHA=6*(0.25-coords(2)*coords(2))
      
      T_USER(1)=0
      T_USER(2)=-1
      T_USER(3)=0
     
      RETURN
      END