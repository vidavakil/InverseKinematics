/**************************************************************************/
//
//                           Disclaimer - Attention
//
//  This software was designed and written by Vida Vakilotojar, as part of
//  her Masters studies in the field of AI & Robotics, at Tehran Polytechnic
//  University. It solves the Inverse Kinematics problem for a robotic arm,
//  and was designed for the MS Windows environment.
//
//  Any portion of this software can be used as-is, only for academic and
//  non-commercial purposes, and provided that this disclaimer is retained
//  exaclty in both source and binary productions or redistributions.
//  If the features of this software are enhanced, or it is ported to
//  other platforms, a version of it should be available to the original
//  author on demand.
//
/**************************************************************************/

#include  <conio.h>
#include  <math.h>
#include  <graphics.h>
#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <time.h>

#define   ITMAX     10                 /* 3*MAXSCREW/CORMAX */
#define   CORMAX    (1.0*M_PI/180.0)   /* degree */
#define   CORMIN    (0.001*M_PI/180.0) /* degree */
#define   PRSMAX    1.0                /* inch   */
#define   PRSMIN    0.01               /* inch   */
#define   ERRMAX    0.005
#define   MAXDISP   1.0                /* inch   */
#define   MAXSCREW  (3.0*M_PI/180.0)   /* degree */
#define   MAXFORCE  2
#define   MAXTRY    20

#define   MAXJOINT       4
#define   DenavitSize    4
#define   SixDegree      6
#define   THREE          3
#define   Zero           0.00000000001

#define   GetAnswer(X)   while(((X=tolower(getch()))!='y')&&(X!='n'));printf("\n")
#define   Sign(X)        ((X>0.0)?1:-1)

typedef   double  DenType[DenavitSize][DenavitSize];
DenType   C, A, Screw, TempMatrix, T_Initial, T_Intermediate,
          T_Final = {{-.64777, -.47174, .598203, -8.9573},
                     {-.14337, -.6957, -.703875, -13.357},
                     {.74822, -.54172, .383020, 12.544},
                     {0, 0, 0, 1}};


typedef   double  CoefType[SixDegree][SixDegree];
CoefType  Coefficients;
double    Constants[SixDegree], Variables[SixDegree], OldVar[SixDegree],
          Accumulator[SixDegree];

CoefType  S_Coeff;
double    S_Const[SixDegree];

enum      JointType { Revolute, Prismatic };
char      Joint[SixDegree] = {Revolute, Revolute, Revolute, Revolute, Revolute, Revolute};

double    FirstAngle[MAXJOINT],
          JointAngle[MAXJOINT];

double    LinkTwist [MAXJOINT] = {M_PI_2, -M_PI_2, M_PI_2, 0},
          LinkLength[MAXJOINT] = {0, 0.375, 0, 5.9},
          LinkOffset[MAXJOINT] = {0, 10, 12.2, 0},
          FirstAngle[MAXJOINT];
char      Robot[30] = "4 DOF";

/*
double    LinkTwist [MAXJOINT] = {-M_PI_2, 0, M_PI_2, -M_PI_2, M_PI_2, 0},
          LinkLength[MAXJOINT] = {0, 431.8/25, -20.32/25, 0, 0, 0},
          LinkOffset[MAXJOINT] = {0, 149.09/25, 0, 433.07/25, 0, 56.25/25},
          JointRange[MAXJOINT][2] ={{-160, 160}, {225, 45}, {-45, 225},
                                    {-110, 170}, {-100, 100}, {-266, 266}};
char      Robot[30] = "PUMA";
*/
double    ScrewAxes[THREE], ScrewPoint[SixDegree], ScrewAngle,
          ScrewTranslation, ScrewArm;

double    P_L, P_M, P_N[3];
double    N_P, O_P, A_P,
          dNx, dNy, dNz, dOx, dOy, dOz, dAx, dAy, dAz, dPx, dPy, dPz;

double    C_Roll, S_Roll, C_Pitch, S_Pitch, C_Yaw, S_Yaw,
          C_Yaw_S_Pitch, S_Yaw_S_Pitch;
double    Y_P_R[3];

double    C_Teta, S_Teta, C_Alfa, S_Alfa;
double    C_Angle, S_Angle, C_Angle_Minus, X_S_Angle, Y_S_Angle, Z_S_Angle;

char      Str[100],
          Vector[4] = {'L', 'M', 'N', 'P'},
          *YPR[3]   = {"YAW", "PITCH", "ROLL"};

int       Count = 0, VibrateCount = 0, Show = 0, TryCount = 0;
double    PositionFactor;


void   CopyMatrix (DenType M1, DenType M2);
void   CrossProduct (DenType T, int N, int L, int M);
void   Identity (DenType Matrix);
void   MulMatrix (DenType M3, DenType M1, DenType M2);
void   Inverse (DenType Inv_Matrix, DenType Matrix);
void   ComputeA (int i);
void   ComputeC (int i);
void   ComputeD (int i);
void   ComputeConstants (DenType TI, DenType TF);
void   Gauss (CoefType Coeff, double *Const, double *Var, int Rank, int Size);
void   ForwardKinematic (DenType Matrix);
void   SolveTeta (DenType TI, DenType TF);
int    ComputeTheta (DenType OldT, DenType NewT);
void   ComputeScrew (DenType Screw);
void   ComputeNewScrew (DenType Screw, int *Divider);
void   ComputeScrewParameters (DenType OldT, DenType NewT, DenType Screw,
                               int *IntermediatePoints);
int    ComputeIntermediatePoint (void);
int    ResolveVibration (int Try);
int    Interpolate (DenType TI, DenType TF);
void   Initialize (DenType Matrix, char *Turn);
void   InitializeGraphics (void);
void   Initialization (void);
void   ComputeT (double P1, double P2, double P3,
                 double Yaw, double Pitch, double Roll,
                 DenType Matrix);
double ComputePosFact (DenType TI, DenType TF);
void   ComputeOrientation (DenType TI);
void   Message (char *S, int Col, int Clear);
void   PrintResult (DenType T, int Column);
int    Converge (DenType TI, DenType TF);

/**************************************************************************/
/* This routine copies the second matrix passesd to it to the first one.  */
/**************************************************************************/
void CopyMatrix (DenType M1, DenType M2)
{
  int i, j;
  for (i = 0; i < DenavitSize; i++)
    for (j = 0; j < DenavitSize; j++)
      M1[i][j] = M2[i][j];
} /* CopyMatrix */


/**************************************************************************/
/*  This routine computes the cross product of the Lth and Mth columns of */
/*  matrix T and places the result into its Nth column; It is called for  */
/*  the formation of the third component of orientation,                  */
/*  eg. N from L and M.                                                   */
/**************************************************************************/
void CrossProduct (DenType T, int N, int L, int M)
{
  T[0][N] = T[1][L] * T[2][M] - T[2][L] * T[1][M];
  T[1][N] = T[2][L] * T[0][M] - T[0][L] * T[2][M];
  T[2][N] = T[0][L] * T[1][M] - T[1][L] * T[0][M];
} /* CrossProduct */


/**************************************************************************/
/*      This routine equates its input matrix as an identity matrix.      */
/**************************************************************************/
void Identity (DenType Matrix)
{
  register int i, j;

  for (i = 0; i < DenavitSize; i++)
    {
      for (j = 0; j < DenavitSize; j++)
        Matrix[i][j] = 0;
      Matrix[i][i] = 1;
    } /* for */
} /* Identity */


/***************************************************************************/
/*  This routine multiplies two denavit_Hartenberg matrices (M1*M2) and    */
/* places the result into M3; it considers the sparseness of D_H matrices. */
/***************************************************************************/
void MulMatrix (DenType M3, DenType M1, DenType M2)
{
  register int col, temp;
  int      row;
  double    Sum;

  for (row = 0; row < DenavitSize - 1 ; row++)
    {
      for (col = 0; col < DenavitSize; col++)
        {
          for (Sum = 0, temp = 0; temp < DenavitSize - 1 ; temp++)
            Sum += M1[row][temp] * M2[temp][col];
          TempMatrix[row][col] = Sum;
        } /* for */
      TempMatrix[row][DenavitSize - 1] += M1[row][DenavitSize - 1];
    } /* for */
  for (row = 0; row < DenavitSize - 1; row++)
    {
      for (col = 0; col < DenavitSize; col++)
        {
          M3[row][col] = TempMatrix[row][col];
        } /* for */
    } /* for */
} /* MulMatrix */


/*****************************************************************************/
/* This routine places the inverse of a D_H matrix (Matrix) into Inv_Matrix. */
/*****************************************************************************/
void Inverse (DenType Inv_Matrix, DenType Matrix)
{
  register int i, j;

  for (i = 0; i < THREE; i++)
    {
      for (j = 0; j < THREE; j++)
        Inv_Matrix[i][j] = Matrix[j][i];
      Inv_Matrix[i][3] = 0;
    } /* for */
  for (i = 0; i < THREE; i++)
    for (j = 0; j < THREE; j++)
      Inv_Matrix[i][3] -= Matrix[j][i] * Matrix[j][3];
} /* Inverse */


/*****************************************************************************/
/* This routine fills the global matrix A with the elements of the ith joint */
/* D_H matrix; this is done by using the fixed joint and link parameters of  */
/* the given robot and the current joint angles.                             */
/* The resultant matrix is used to compute Ci where :                        */
/*      Ci = A0*A1*...*Ai = Ci-1*Ai        i = 0,1,...5  for a 6_dof robot   */
/*****************************************************************************/
void ComputeA (int i)
{

  C_Teta = cos (JointAngle[i]);
  S_Teta = sin (JointAngle[i]);
  C_Alfa = cos (LinkTwist[i]);   /* fix value */
  S_Alfa = sin (LinkTwist[i]);   /* fix value */
  A[0][0] = C_Teta;
  A[0][1] = - S_Teta * C_Alfa;
  A[0][2] = S_Alfa * S_Teta;
  A[1][0] = S_Teta;
  A[1][1] = C_Teta * C_Alfa;
  A[1][2] = - C_Teta * S_Alfa;
  A[2][0] = 0;
  A[2][1] = S_Alfa;
  A[2][2] = C_Alfa;
  A[0][3] = LinkLength[i] * C_Teta;
  A[1][3] = LinkLength[i] * S_Teta;
  A[2][3] = LinkOffset[i];
} /* ComputeA */


/***************************************************************************/
/*    This routine computes Ci as described previously :   Ci = Ci-1*Ai    */
/*    and places the result into the global matrix C.                      */
/*    if i = 0 then A0 is copied to C :    C0 = A0                         */
/***************************************************************************/
void ComputeC (int i)
{
  register int j, k;

  ComputeA (i);
  if (i == 0)
    CopyMatrix (C, A);
  else
    MulMatrix (C, C, A);
} /* ComputeC */


/******************************************************************************/
/* This routine computes the necessary elements of the ith column of the      */
/* Jacobian matrix. Actually it computes the necessary elements of the Di     */
/* matrix where    Di = Ci-1 * Qi * (Ci-1)Inv                                 */
/* and depending on the type of the ith joint we have                         */
/*                                                                            */
/*          | 0  -1   0   0 |                                                 */
/*          | 1   0   0   0 |                                                 */
/*     Qi = | 0   0   0   0 |  the ith joint is revolute                      */
/*          | 0   0   0   0 |                                                 */
/*                                                                            */
/*          | 0   0   0   0 |                                                 */
/*          | 0   0   0   0 |                                                 */
/*     Qi = | 0   0   0   1 |  the ith joint is prismatic                     */
/*          | 0   0   0   0 |                                                 */
/*                                                                            */
/* where                                                                      */
/*     Qi = Bi*Ai_Inv                                                         */
/* and                                                                        */
/*     Bi = partial derivative of Ai with respect to the ith joint angle      */
/*          computed at current angles                                        */
/*                                                                            */
/* therefore Di takes a simple form                                           */
/*                                                                            */
/*          |  0     -N3      N2      (P cross N)1 |                          */
/*          |  N3     0      -N1      (P cross N)2 |                          */
/*     Di = | -N2     N1      0       (P cross N)3 | the ith joint is revolute*/
/*          |  0      0       0             0      |                          */
/*                                                                            */
/*          |  0      0       0       N1 |                                    */
/*          |  0      0       0       N2 |                                    */
/*     Di = |  0      0       0       N3 | the ith joint is prismatic         */
/*          |  0      0       0       0  |                                    */
/*                                                                            */
/* where we have                                                              */
/*                                                                            */
/*          | L1   M1   N1   P1 |                                             */
/*          | L2   M2   N2   P2 |                                             */
/*   Ci-1 = | L3   M3   N3   P3 |                                             */
/*          | 0    0    0    1  |                                             */
/*                                                                            */
/* Notice : D0 = Q0                                                           */
/*                                                                            */
/* The elements of the D={dij} matrices that have been used to form the       */
/* Jacobian are : d14, d24, d34, d12, d13 and d23 respectively.               */
/* For lower DOF robots the remaining columns of Jacobian are not filled, but */
/* all of the six equations relating to position and orientation are used to  */
/* solve for the joint changes; Therefore for a four DOF robot, unless there  */
/* is no solution at all, all of the variables will be computed by the first  */
/* four suitable equations out of the six ones.                               */
/******************************************************************************/
void ComputeD (int i)
{
  register int j;

  if (i == 0)
    {
      if (Joint[0] == Revolute)
        {
          for (j = 0; j < MAXJOINT; j++)
            Coefficients[j][0] = 0;
          Coefficients[3][0] = -1;
        } /* if */
      else
        {
          for (j = 0; j < MAXJOINT; j++)
            Coefficients[j][0] = 0;
          Coefficients[2][0] = 1;
        } /* else */
    } /* if */
  else
    {
      ComputeC (i - 1);
      if (Joint[i] == Revolute)
        {
          P_N[0] = C[1][3] * C[2][2] - C[2][3] * C[1][2];
          P_N[1] = C[2][3] * C[0][2] - C[0][3] * C[2][2];
          P_N[2] = C[0][3] * C[1][2] - C[1][3] * C[0][2];
          Coefficients[0][i] = P_N[0];     /* (P cross N)(1) */
          Coefficients[1][i] = P_N[1];     /* (P cross N)(2) */
          Coefficients[2][i] = P_N[2];     /* (P cross N)(3) */
          Coefficients[3][i] = - C[2][2];   /* -N3 */
          Coefficients[4][i] = C[1][2];     /* +N2 */
          Coefficients[5][i] = - C[0][2];   /* -N1 */
        } /* if */
      else
        {
          for (j = 0; j < 3; j++)
            Coefficients[j][i] = C[j][2];
          for (j = 3; j < MAXJOINT; j++)
            Coefficients[j][i] = 0;
        } /* else */
    } /* else */
} /* ComputeD */


/***************************************************************************/
/* This routine computes the necessary elements of   (TF-TI)*TI_Inv        */
/* which are needed to solve the following equation for the changes of     */
/* joint angles :                                                          */
/*                                                                         */
/*                             __ 6                                        */
/*            (TF-TI)*TI_Inv = \ Di*DeltaTeta                              */
/*                             /_i=1                                       */
/* or                                                                      */
/*            T    = J   * Teta                                            */
/*            =6*1   =6*6   =6*1                                           */
/* The six elements of (TF-TI)*TI_Inv are selected as for Di's.            */
/* N_P , O_P and A_P are inner products of the associated columns of TF.   */
/*                                                                         */
/***************************************************************************/
void ComputeConstants (DenType TI, DenType TF)
{
  register int j;

  N_P = O_P = A_P = 0;
  for (j = 0; j < 3; j++)
    {
      N_P += TF[j][0] * TF[j][3];
      O_P += TF[j][1] * TF[j][3];
      A_P += TF[j][2] * TF[j][3];
    } /* for */
  dNx = TF[0][0] - TI[0][0];
  dNy = TF[1][0] - TI[1][0];
  dNz = TF[2][0] - TI[2][0];
  dOx = TF[0][1] - TI[0][1];
  dOy = TF[1][1] - TI[1][1];
  dOz = TF[2][1] - TI[2][1];
  dAx = TF[0][2] - TI[0][2];
  dAy = TF[1][2] - TI[1][2];
  dAz = TF[2][2] - TI[2][2];
  dPx = TF[0][3] - TI[0][3];
  dPy = TF[1][3] - TI[1][3];
  dPz = TF[2][3] - TI[2][3];
  Constants [0] = - dNx * N_P - dOx * O_P - dAx * A_P + dPx;
  Constants [1] = - dNy * N_P - dOy * O_P - dAy * A_P + dPy;
  Constants [2] = - dNz * N_P - dOz * O_P - dAz * A_P + dPz;
  Constants [3] = dNx * TI[1][0] + dOx * TI[1][1] + dAx * TI[1][2];
  Constants [4] = dNx * TI[2][0] + dOx * TI[2][1] + dAx * TI[2][2];
  Constants [5] = dNy * TI[2][0] + dOy * TI[2][1] + dAy * TI[2][2];
} /* ComputeConstants */


/****************************************************************************/
/* This routine solves a system of linear equations using the method of     */
/* Gauss. It is called to solve for the S and D parameters of a             */
/* Screw_Displacement matrix, and also to solve the Jacobian equation for   */
/* the change of angles. Here Size is the number of equations and Rank is   */
/* the number of Variables to be solved. If Rank is smaller than Size,      */
/* during the pivoting of the coefficients matrix, the set of proper        */
/* eqauations by means of which the system can be solved is extracted and   */
/* put on the first Rank rows of the coefficients matrix; then the Rank by  */
/* Rank system of equations can be solved easily for the demanded Variables.*/
/* If a column of the Coeff matrix (matrix of coefficients) is zero then    */
/* the corresponding variable is irrelenvant and the system has no special  */
/* solution; the program exits with an error message.                       */
/* If two variables of the system are dependent, it is detected, and one of */
/* them is set to zero to solve the system; then using the value gained for */
/* the other variable, the previous values of the joint angle changes, and  */
/* the coefficients of the variables in the system of equations,            */
/* both of the variables are computed such that the ratio and the sign of   */
/* the changes remains the same. In this way degeneracies can be resolved.  */
/* This condition does not happen when solving for Screw parameters.        */
/* For each joint the previous changes are retained in OldVar matrix,       */
/* and the current changes are put into Variables matrix.                   */
/****************************************************************************/
void Gauss (CoefType Coeff, double *Const, double *Var, int Rank, int Size)
{
  double Temp, C1, C2;
  int   i, j, k, Couple2 = -1, Couple1 = -1;

  for (k = 0; k < Rank; k++)
    {
      for (i = k + 1, j = k, Temp = fabs (Coeff[k][k]); i < Size; i++)
        if (fabs (Coeff[i][k]) > Temp)
          {
            Temp = fabs (Coeff[i][k]);
            j = i;
          } /* if */
      if (Temp < Zero)
        {
          for (i = k - 1; i >= 0; i--)
            if (fabs (Coeff[i][k]) > Zero)
              {
                Couple1 = i;
                break;
              } /* if */
          if (Couple1 == -1)
            {
              Message ("No Special Solution ", 420, 1);
              getch ();
              closegraph ();
              exit (-1);
            } /* if */
          Couple2 = k;
          sprintf (Str, "Coupling of Joints %d & %d", Couple1+1, Couple2+1);
          Message (Str, 420, 1);
          for (i = Couple2+1; i < Rank; i++)
            Coeff[Couple2][i] = 0;
          Coeff[Couple2][Couple2] = 1.0;
          Const[Couple2] = 0;
          if (Couple1 != -1)
            {
              C1 = Coeff[Couple1][Couple1];
              C2 = Coeff[Couple1][Couple2];
            } /* if */
          } /* if */
      else if (j != k)
        {
          for (i = k; i < Rank; i++)
            {
              Temp = Coeff[k][i];
              Coeff[k][i] = Coeff[j][i];
              Coeff[j][i] = Temp;
            } /* for */
          Temp = Const[k];
          Const[k] = Const[j];
          Const[j] = Temp;
        } /* if */
      for (i = k + 1; i < Size; i++)
        {
          Temp = (double) Coeff[i][k] / (double) Coeff[k][k];
          for (j = k; j < Size; j++)
            Coeff[i][j] -= Coeff[k][j] * Temp;
          Const[i] -= Const[k] * Temp;
        } /* for */
    } /* for */
  for (i = Rank - 1; i >= 0; i--)
    {
      Temp = 0;
      for (Temp = 0, j = i + 1; j < Size; j++)
        Temp += Var[j] * Coeff[i][j];
      Var[i] = (double) (Const[i] - Temp) / (double) Coeff[i][i];
    } /* for */

  if ((Couple2 != -1) && (Couple1 != -1))
    {
      if ((fabs (OldVar[Couple2]) > Zero) && (fabs (OldVar[Couple1]) > Zero))
        {
          Variables[Couple2] = Variables[Couple1] * C1 * Sign (OldVar[Couple2]) /
                           (C1 * Sign (OldVar[Couple1]) / OldVar[Couple2] + C2);
          Variables[Couple1] = (Variables[Couple1] * C1 - C2 * Variables[Couple2])/C1;
        } /* if */
      else
        if (fabs (OldVar[Couple2]) > Zero)
          {
            Variables[Couple2] = Variables[Couple1] * C1 / C2;
            Variables[Couple1] = 0;
          } /* else */
    } /* if */
} /* Gauss */


/****************************************************************************/
/* This routine solves for Forward_Kinematic using the current joint angles */
/* by multiplying the D_H matrices of joints, and placing the result in     */
/* Matrix.                                                                  */
/****************************************************************************/
void ForwardKinematic (DenType Matrix)
{
  register int i;

  Identity (Matrix);
  for (i = 0; i < MAXJOINT; i++)
    {
      ComputeA (i);
      MulMatrix (Matrix, Matrix, A);
    } /* for */
} /* ForwardKinematic */


/****************************************************************************/
/* This routine forms the Jacobian Matrix and equation and then calls Gauss */
/* to solve for the joint angle changes. Althouh the given robot might      */
/* have less than 6 DOF, but the routine Gauss is called with SixDegree     */
/* instead of MAXJOINT to let the Joint angle changes be computed with more */
/* accuracy and less risk.                                                  */
/****************************************************************************/
void SolveTeta (DenType TI, DenType TF)
{
  register int i;

  for (i = 0; i < MAXJOINT; i++)
    {
      ComputeD (i);
      OldVar[i] = Variables[i];
    } /* for */
  ComputeConstants (TI, TF);
  Gauss (Coefficients, Constants, Variables, MAXJOINT, SixDegree);
} /* SolveTeta */

/*****************************************************************************/
/* This routine first calls SolveTeta to compute the changes of the joint    */
/* angles. Then it finds the max ratio between a change and the max allowable*/
/* change (for revolute joints using CORMAX and for prismatic joints using   */
/* PRSMAX). If this ratio is greater than one then all the changes are       */
/* divided by this ratio. Then if a change is less than the minimum          */
/* allowable change, it is rounded either to 0 or to this minimum change.    */
/* Then this changes are applied to joint angles, while retaining them in    */
/* the range -2*PI to 2*PI. At the same time a range check of the angles can */
/* be done. If the accumulated changes of an angle reaches CORMAX or PRSMAX  */
/* then the new position and orientation will be printed on the screen by    */
/* calling PrintResult. If none of the joint angles have changed this is     */
/* reported by returning a 0, else a 1 is returned; therefore 0 denotes      */
/* that no more changes can be reached through the current iteration.        */
/*****************************************************************************/
int ComputeTheta (DenType TI, DenType TF)
{
  register int i, j;
  int      Div = 1;
  double   Ratio, TempRatio;

  Count++;
  SolveTeta (TI, TF);
  Ratio = fabs (Variables[0] / (double) ((Joint[0] == Revolute)
                                                     ? CORMAX : PRSMAX));
  for (i = 1; i < MAXJOINT; i++)
    if ((TempRatio = fabs (Variables[i] / (double) ((Joint[i] == Revolute)
                                         ? CORMAX : PRSMAX))) > Ratio)
      Ratio = TempRatio;
  if (Ratio > 1.0)
    for (i = 0; i < MAXJOINT; i++)
      Variables[i] /= (double) Ratio;
  for (i = 0; i < MAXJOINT; i++)
    {
      Ratio = (Joint[i] == Revolute) ? CORMIN : PRSMIN; /* not a real ratio */
      if (fabs (Variables[i]) < Ratio)
        {
          if (fabs (Variables[i]) > Ratio / 2.0)
            Variables[i] = Ratio * ((Variables[i] > 0.0) ? 1 : -1);
          else
            Variables[i] = 0;
        } /* if */
    } /* for */
  for (i = 0, Ratio = 0; i < MAXJOINT; i++)
    Ratio += fabs (Variables[i]);
  for (i = 0; i < MAXJOINT; i++)
    {
      JointAngle[i] += Variables[i];
      Accumulator[i] += Variables[i];
      if (fabs (Accumulator[i]) > CORMAX)
        Show = 1;
      if ((Joint[i] == Revolute) && (fabs (JointAngle[i]) > (2.0 * M_PI)))
        {
          Div = JointAngle[i] / (2.0 * M_PI);
          JointAngle[i] -= Div * 2.0 * M_PI;
        } /* if */
/*
      if (JointAngle[i] < JointRange[i][0] * M_PI / 180)
        {
          JointAngle[i] = JointRange[i][0] * M_PI / 180;
          Variables[i] = 0;
        }
      if (JointAngle[i] > JointRange[i][1] * M_PI / 180)
        {
          JointAngle[i] = JointRange[i][1] * M_PI / 180;
          Variables[i] = 0;
        }
*/
    } /* for */
  ForwardKinematic (TI);
  if (Show)
    {
      PrintResult (TI, 20);
      for (i = 0; i < MAXJOINT; i++)
        Accumulator[i] = 0;
      Show = 0;
    } /* if */
  if (Ratio != 0)
    return (0);
  else
    return (1);
} /* ComputeTheta */


/******************************************************************************/
/* This routine computes a Screw_Displacement matrix into Screw by using      */
/* the global variables : ScrewAngle, ScrewTranslation, ScrewPoint, and       */
/* ScrewAxis. After the Screw_Displacement required to move the tool frame    */
/* from the current position and orientation to the final point is computed   */
/* the ScrewAngle and ScrewTranslation are broken to smaller ones; noticing   */
/* that the length of a helix arrond an axis with a radius of ScrewArm and    */
/* an angle of Screwangle and hight of ScrewTranslation is :                  */
/*   ((ScrewTranslation^2) + (ScrewArm*ScrewAngle)^2) ^(1/2)                  */
/* This is compared to MAXDISP, and ScrewAngle is compared to MAXSCREW.       */
/* Then by calling this routine a new Screw_Displacement matrix is generated  */
/* arround the same axis with the new smaller ScrewAngle and ScrewTranslation */
/******************************************************************************/
void ComputeScrew (DenType Screw)
{
  register int i, j;

  C_Angle = cos (ScrewAngle);
  C_Angle_Minus = 1.0 - C_Angle;
  S_Angle = sin (ScrewAngle);
  for (i = 0; i < THREE; i++)
    {
      for (j = 0; j < THREE; j++)
        Screw[i][j] = ScrewAxes[i] * ScrewAxes[j] * C_Angle_Minus;
      Screw[i][i] += C_Angle;
    } /* for */
  X_S_Angle = ScrewAxes[0] * S_Angle;
  Y_S_Angle = ScrewAxes[1] * S_Angle;
  Z_S_Angle = ScrewAxes[2] * S_Angle;
  Screw[0][1] -= Z_S_Angle;
  Screw[1][2] -= X_S_Angle;
  Screw[2][0] -= Y_S_Angle;
  Screw[0][2] += Y_S_Angle;
  Screw[1][0] += Z_S_Angle;
  Screw[2][1] += X_S_Angle;
  for (i = 0; i < THREE; i++)
    {
      Screw[i][3] = ScrewTranslation * ScrewAxes[i] + ScrewPoint[i];
      for (j = 0; j < THREE; j++)
        Screw[i][3] -= ScrewPoint[j] * Screw[i][j];
    } /* for */
  for (i = 0; i < THREE; i++)
    Screw[3][i] = 0;
  Screw[3][3] = 1;
} /* ComputeScrew */


/******************************************************************************/
/* As described in the previous comment, by considering the Screw parameters, */
/* and MAXDISP and MAXSCREW, new ScrewAngle and ScrewTranslation are computed */
/* and then ComputeScrew is called to generate the corresponding Screw matrix */
/* Here Divider is an integer that initially is passed to this routine to be  */
/* used for the division of ScrewAngle and ScrewTranslation. It is computed   */
/* by just considering the euclidean distance of initial and final positions. */
/* If it is less than required it will be changed properly by this routine.   */
/******************************************************************************/
void ComputeNewScrew (DenType Screw, int *Divider)
{
  double Ratio1, Ratio;

  Ratio = sqrt (pow (ScrewTranslation, 2) + pow (ScrewArm * ScrewAngle, 2)) / MAXDISP;
  if ((Ratio1 = fabs (ScrewAngle) / MAXSCREW) > Ratio)
    Ratio = Ratio1;
  if (*Divider > (int) Ratio)
    Ratio = (double) *Divider;
  else
    *Divider = (int) Ratio;
  ScrewTranslation /= (double) Ratio;
  ScrewAngle /= (double) Ratio;
  ComputeScrew (Screw);
} /* ComputeNewScrew */


/******************************************************************************/
/* This routine takes the initial(OldT) and final(NewT) frames and computes   */
/* the Screw_Displacement matrix required to move the initial frame to the    */
/* final one as follows   :  Screw = Newt*OldT_Inv                            */
/* Then it computes the following parameters of this Screw_Displacement :     */
/* ScrewAxis, ScrewAngle, ScrewPoint, ScrewTranslation and  ScrewArm.         */
/* Screwpoint is the intersection point of the ScrewAxis with a               */
/* perpendicular plane that passes through the first point (OldT).            */
/* ScrewArm is the radius of rotation about the ScerwAxis. ScrewAxis is the   */
/* vector parallel to axis of screw. ScrewTranslation is the translation along*/
/* the ScrewAxis. If the ScrewAngle is too small (its cos is near to 1),      */
/* it is set to zero and Screwpoint is set to the first point (OldT),         */
/* ScrewAxes is computed by dividing the fourth column elements of Screw by   */
/* ScrewTranslation which is the length of that column, and if                */
/* ScrewTranslation is too small ScrewAngle is set to be parallel to the Z    */
/* axis.                                                                      */
/* If ScrewAngle is not too small ScrewAxis is computed by the procedures     */
/* described in PAUL's book.                                                  */
/* Then a linear system of equations for finding ScrewPoint and               */
/* ScrewTranslation is generated by considering the equations of the fourth   */
/* column of the Screw_Displacement matrix and another equation that selects  */
/* the ScrewPoint to be the point of intersection of the ScrewAxis with a     */
/* perpendicular plane passing through the first point of Screw(OldT).        */
/* Then Gauss is called to solve this system of equations. And Now ScrewArm   */
/* can be computed, and ComputeNewScew can be called to generate an           */
/* intermediate Screw to move the initial frame to an intermediate frame,     */
/* closer to the final frame.                                                 */
/******************************************************************************/
void ComputeScrewParameters (DenType OldT, DenType NewT, DenType Screw,
                             int *IntermediatePoints)
{
  double S_Angle, Vers, Max;
  register int i, j;

  Inverse (Screw, OldT);
  MulMatrix (Screw, NewT, Screw);
  S_Angle = (Screw[0][0] + Screw[1][1] + Screw[2][2] - 1.0) / 2.0;
  if (fabs (S_Angle) > 0.9999)
    {
      ScrewAngle = 0;
      for (i = 0, ScrewTranslation = 0; i < THREE; i++)
        {
          ScrewPoint[i] = OldT[i][3];
          ScrewTranslation += pow (Screw[i][3], 2);
        } /* for */
      ScrewTranslation = sqrt (ScrewTranslation);
      if (ScrewTranslation > 0.00000001) /* 0.000001 */
        for (i = 0; i < THREE; i++)
          ScrewAxes[i] = Screw[i][3] / ScrewTranslation;
      else
        {
          ScrewAxes[0] = ScrewAxes[1] = 0;
          ScrewAxes[2] = 1;
        } /* else */
    } /* for */
  else
    {
/*    ScrewAngle = (double) acos ((double) S_Angle);*/
      ScrewAngle = (double) atan2 (sqrt (pow (Screw[2][1] - Screw[1][2], 2) +
                                         pow (Screw[2][0] - Screw[0][2], 2) +
                                         pow (Screw[1][0] - Screw[0][1], 2)),
                                        (S_Angle * 2.0));

      if (ScrewAngle > M_PI_2)
        {
          Vers = 1 - S_Angle;
          ScrewAxes[0] = Sign (Screw[2][1] - Screw[1][2]) *
                         sqrt ((Screw[0][0] - S_Angle) / Vers);
          ScrewAxes[1] = Sign (Screw[0][2] - Screw[2][0]) *
                         sqrt ((Screw[1][1] - S_Angle) / Vers);
          ScrewAxes[2] = Sign (Screw[1][0] - Screw[0][1]) *
                         sqrt ((Screw[2][2] - S_Angle) / Vers);
          for (i = j = 0, Max = Screw[0][0]; i < 3; i++)
            if (Screw[i][i] > Max)
              {
                Max = Screw[i][i];
                j = i;
              } /* if */
          Max = 2 * ScrewAxes[j] * Vers;
          switch (j)
            {
              case 0 : ScrewAxes[1] = (Screw[1][0] + Screw[0][1]) / Max;
                       ScrewAxes[2] = (Screw[2][0] + Screw[0][2]) / Max;
                       break;
              case 1 : ScrewAxes[0] = (Screw[1][0] + Screw[0][1]) / Max;
                       ScrewAxes[2] = (Screw[1][2] + Screw[2][1]) / Max;
                       break;
              case 2 : ScrewAxes[0] = (Screw[0][2] + Screw[2][0]) / Max;
                       ScrewAxes[1] = (Screw[1][2] + Screw[2][1]) / Max;
                       break;
            } /* switch */
        } /* if */
      else
        {
          S_Angle = 2 * sin (ScrewAngle);
          ScrewAxes[0] = (Screw[2][1] - Screw[1][2]) / S_Angle;
          ScrewAxes[1] = (Screw[0][2] - Screw[2][0]) / S_Angle;
          ScrewAxes[2] = (Screw[1][0] - Screw[0][1]) / S_Angle;
        } /* else */
      for (i = 0, S_Const[3] = 0; i < THREE; i++)
        {
          S_Coeff[i][0] = - Screw[i][0];
          S_Coeff[i][1] = - Screw[i][1];
          S_Coeff[i][2] = - Screw[i][2];
          S_Coeff[i][3] = S_Coeff[3][i] = ScrewAxes[i];
          S_Const[i] = Screw[i][3];
          S_Const[3] += OldT[i][3] * ScrewAxes[i];
        } /* for */
      S_Coeff[0][0] += 1;
      S_Coeff[1][1] += 1;
      S_Coeff[2][2] += 1;
      S_Coeff[3][3] = 0;
      Gauss (S_Coeff, S_Const, ScrewPoint, 4, 4);
      ScrewTranslation = ScrewPoint[3];
    } /* else */
  for (i = 0, ScrewArm = 0; i < 3; i++)
    ScrewArm += pow ((OldT[i][3] - ScrewPoint[i]), 2);
  ScrewArm = sqrt (ScrewArm);
  ComputeNewScrew (Screw, IntermediatePoints);
} /* ComputeScrewParameters */


/******************************************************************************/
/* This routine first computes number of intermediate points required to      */
/* move the initial frame to the final one by just considering the euclidean  */
/* distance of the two frames. Then it calls ComputeScrewParameters with this */
/* number, which is changed properly by that routine and a new Screw is       */
/* generated to move the initial frame to an intermediate one closer to the   */
/* final. If IntermediatePoints is 0 then convergence is gained, else         */
/* the computed Screw matrix is multiplied to the initial frame to generate   */
/* the intermediate frame. Finally IntermediatePoints is returned to the      */
/* caller                                                                     */
/******************************************************************************/
int ComputeIntermediatePoint (void)
{
  register int i;
  int      IntermediatePoints;
  double    Distance;

  for (i = 0, Distance = 0; i < 3; i++)
    Distance += pow ((T_Initial[i][3] - T_Final[i][3]), 2);
  IntermediatePoints = sqrt (Distance) / MAXDISP;
  ComputeScrewParameters (T_Initial, T_Final, Screw, &IntermediatePoints);
  if (!IntermediatePoints)
    return (0);
  CopyMatrix (T_Intermediate, T_Initial);
  MulMatrix (T_Intermediate, Screw, T_Intermediate);
  return (1);
} /* ComputeIntermediatePoint */


/**************************************************************************/
/* This routine considers the result of the last call to Interpolate.     */
/* If it was successful it just copies current angles in FirstAngle and   */
/* sets TryCount to 0. Else if TryCount has reached MAXTRY no more        */
/* convergence can be gained at all and a 1 is returned to tell the caller*/
/* to exit. Otherwise it looks to see which joint has the greatest change */
/* since the last time by computing the difference between JointAngle and */
/* Firstangle. Then it forces that joint to move for MAXFORCE times, each */
/* time at an angle of CORMAX, in the direction of its last changes.      */
/* Each time the resultant frame is printed. In this way vibrations       */
/* arround a point can be got rid of. Finally FirstAngle is updated to    */
/* current joint angles.                                                  */
/**************************************************************************/
int ResolveVibration (int Try)
{
  register int i, j;
  double    MaxChange, Change;
  int       Max;

  if (!Try)
    {
      for (i = 0; i < MAXJOINT; i++)
        FirstAngle[i] = JointAngle[i];
      TryCount = 0;
      return (0);
    } /* if */

  if (TryCount++ == MAXTRY)
    return (1);

  VibrateCount++;
  MaxChange = fabs (JointAngle[0] - FirstAngle[0]);
  for (i = 0, Max = 0; i < MAXJOINT; i++)
    {
      if ((Change = fabs (JointAngle[i] - FirstAngle[i])) > MaxChange)
        {
          MaxChange = Change;
          Max = i;
        } /* if */
    } /* for */
  j = Sign (JointAngle[Max] - FirstAngle[Max]);
  for (i = 0; i < MAXFORCE; i++)
    {
      JointAngle[Max] += j * CORMAX;
/*
      if (JointAngle[Max] < JointRange[Max][0]*M_PI/180)
        {
          JointAngle[Max] = JointRange[Max][0]*M_PI/180;
          i = MAXFORCE;
        }
      if (JointAngle[Max] > JointRange[Max][1]*M_PI/180)
        {
          JointAngle[Max] = JointRange[Max][1]*M_PI/180;
          i = MAXFORCE;
        }
*/
      ForwardKinematic (T_Initial);
      PrintResult (T_Initial, 20);
      sprintf (Str, "Vibration Resolution of Joint%d", Max+1);
      Message (Str, 420, 1);
    } /* for */
  for (i = 0; i < MAXJOINT; i++)
    FirstAngle[i] = JointAngle[i];
  return (0);
} /* ResolveVibration */


/******************************************************************************/
/* This routine is called to make the initial frame move to the next          */
/* intermediate frame. For ITMAX maximum number of iterations, or until       */
/* convergence of the initial and intermediate frames, ComputeTheta is called.*/
/* In the case of no more convergence, full convergence, or getting out of    */
/* iteration proper messages are printed, and only in the latter case a one   */
/* is returned; if this case happens for MAXTRY consecutive times,            */
/* NO_MORE_CONVERGENCE will be concluded.                                     */
/* The total number of iterations will also be printed.                       */
/******************************************************************************/
int Interpolate (DenType TI, DenType TF)
{
  register int i, ReturnValue;

  PrintResult (TF, 120);
  for (i = 0; (i < ITMAX) && !Converge (TI, TF); i++)
    if (ComputeTheta (TI, TF))
      {
        Message ("No more convergence", 420, 1);
        break;
      } /* if */
  PrintResult (TI, 20);
  if (Converge (TI, TF))
    {
      sprintf (Str, "Done in %d steps", i+1);
      Message (Str, 420, 1);
      return (0);
    } /* else */
  if (i == ITMAX)
    {
      Message ("Out of Iteration", 420, 1);
      return (1);
    } /* if */
  else
    {
      sprintf (Str, "No More Convergence %d steps", i+1);
      Message (Str, 420, 1);
      return (0);
    } /* else */

} /* Interpolate */


/**************************************************************************/
/* This routine is called to generate the initial or finall frames by     */
/* asking questions about the position and orientation of that frames.    */
/* Angles can be given in both radian and degree.                         */
/* Frames can be defined by giving joint angles.                          */
/* In defining the orientation both (L,M,N) and (Yaw_Pitch_Roll) can be   */
/* given.                                                                 */
/**************************************************************************/
void Initialize (DenType Matrix, char *Turn)
{
  register int i;
  char     IfJoint, IfRadian, IfL, IfM;
  double    Length;

  clrscr ();
  printf ("\nThrough the following steps enter %s point.\n", Turn);
  printf ("Would you like to enter %s point parameters in joint coordinate ? (Y/N)", Turn);
  GetAnswer (IfJoint);
  if (IfJoint == 'y')
    {
      printf ("Would you enter joint angles in radian ? (Y/N)");
      GetAnswer (IfRadian);
      printf ("Enter %s point joint parameters.\n", Turn);
      for (i = 0; i < MAXJOINT; i++)
        {
          printf ("Joint[%d] = ", i+1);
          scanf ("%lf", &JointAngle[i]);
          if ((IfRadian == 'n') && (Joint[i] == Revolute))
            JointAngle[i] *= M_PI / 180.0;
        } /* for */
      ForwardKinematic (Matrix);
    } /* if */
  else
    {
      printf ("Enter %s point position parameters.\n", Turn);
      for (i = 0; i < THREE; i++)
        {
          printf ("P[%d] = ", i+1);
          scanf ("%lf", &Matrix[i][3]);
        } /* for */
      printf ("Would you enter Y_P_R parameters ? (Y/N)");
      GetAnswer (IfL);
      if (IfL == 'y')
        {
          printf ("Would you enter Y_P_R angles in radian ? (Y/N)");
          GetAnswer (IfRadian);
          printf ("Enter %s Orientation parameters.\n", Turn);
          for (i = 0; i < THREE; i++)
            {
              printf ("%s = ", YPR[i]);
              scanf ("%lf", &Y_P_R[i]);
              if (IfRadian == 'n')
                Y_P_R[i] *= M_PI / 180.0;
            } /* for */
          ComputeT (Matrix[0][3], Matrix[1][3], Matrix[2][3],
                     Y_P_R[0], Y_P_R[1], Y_P_R[2], Matrix);
        } /* if */
      else
        {
          printf ("Would you enter 'l' parameters ? (Y/N)");
          GetAnswer (IfL);
          if (IfL == 'y')
            {
              for (Length = 0, i = 0; i < THREE; i++)
                {
                  printf ("L[%d] = ", i+1);
                  scanf ("%lf", &Matrix[i][0]);
                  Length += pow (Matrix[i][0], 2);
                } /* for */
              Length = sqrt (Length);
              for (i = 0; i < THREE; i++)
                Matrix[i][0] /= Length;
            } /* if */
          printf ("Would you enter 'm' parameters ? (Y/N)");
          GetAnswer (IfM);
          if (IfM == 'y')
            {
              for (Length = 0, i = 0; i < THREE; i++)
                {
                  printf ("M[%d] = ", i+1);
                  scanf ("%lf", &Matrix[i][1]);
                  Length += pow (Matrix[i][1], 2);
                } /* for */
              Length = sqrt (Length);
              for (i = 0; i < THREE; i++)
                Matrix[i][1] /= Length;
            } /* if */
          if ((IfL == 'n') && (IfM == 'n'))
            {
              printf ("Two of the three vectors l, m, n must be given.");
              exit (-1);
            } /* if */
          if ((IfL == 'y') && (IfM == 'y'))
            CrossProduct (Matrix, 2, 0, 1);
          else
            {
              for (Length = 0, i = 0; i < THREE; i++)
                {
                  printf ("N[%d] = ", i+1);
                  scanf ("%lf ", &Matrix[i][2]);
                  Length += pow (Matrix[i][2], 2);
                } /* for */
              Length = sqrt (Length);
              for (i = 0; i < THREE; i++)
                Matrix[i][2] /= Length;
              if (IfL != 'y')
                CrossProduct (Matrix, 0, 1, 2);
              else
                CrossProduct (Matrix, 1, 2, 0);
              for (i = 0; i < THREE; i++)
                Matrix[3][i] = 0;
            } /* else */
          Matrix[3][3] = 1;
        } /* else */
    } /* else */
  getch();
} /* Initialize */


/**************************************************************************/
/* This routine initializes the graphic screen and draws its windows      */
/**************************************************************************/
void InitializeGraphics (void)
{
   int gdriver = DETECT, gmode, errorcode;
   initgraph(&gdriver, &gmode, "");
   errorcode = graphresult();
   if (errorcode != grOk)  /* an error occurred */
     {
        printf("Graphics error: %s\n", grapherrormsg(errorcode));
        printf("Press any key to halt:");
        getch();
        exit(1);             /* return with error code */
     } /* if */
  setcolor (LIGHTRED);
  rectangle (0, 0, 630, 98);
  outtextxy (310, 40 , " Current Position & Orientation & Y_P_R");
  setcolor (LIGHTBLUE);
  rectangle (0, 100, 630, 198);
  outtextxy (310, 140, " Next Position & Orientation & Y_P_R");
  setcolor (LIGHTCYAN);
  rectangle (0, 200, 630, 298);
  outtextxy (310, 240, " Final Position & Orientation & Y_P_R");
  setcolor (LIGHTGREEN);
  rectangle (0, 300, 630, 398);
  sprintf (Str, " Current Joint Angles of %s", Robot);
  outtextxy (150, 320, Str);
  setcolor (YELLOW);
  rectangle (0, 400, 630, 479);
  outtextxy (260, 420, " Messages");
  setcolor (WHITE);
  setfillstyle (SOLID_FILL, BLACK);
} /* InitializeGraphics */


/**************************************************************************/
/* This routine initializes the initial and final frames, initiaizes      */
/* the screen, prints the final frame, and initializes FirstAngle.        */
/**************************************************************************/
void Initialization (void)
{
  register int i;

  Initialize (T_Final, "Final");
  Initialize (T_Initial, "initial");

  InitializeGraphics ();
  PrintResult (T_Final, 220);

  for (i = 0; i < MAXJOINT; i++)
    FirstAngle[i] = JointAngle[i];


} /* Initialization */


/******************************************************************************/
/* This routine computes a frame matrix given the position and Yaw_Pitch_Roll */
/* parameters. It is called by Initialize.                                    */
/******************************************************************************/
void ComputeT (double P1, double P2, double P3,
                     double Yaw, double Pitch, double Roll,
                     DenType Matrix)
{
  C_Roll = cos (Roll);
  S_Roll = sin (Roll);
  C_Pitch = cos (Pitch);
  S_Pitch = sin (Pitch);
  C_Yaw = cos (Yaw);
  S_Yaw = sin (Yaw);
  C_Yaw_S_Pitch = C_Yaw * S_Pitch;
  S_Yaw_S_Pitch = S_Yaw * S_Pitch;

  Matrix[0][0] = C_Yaw * C_Pitch;
  Matrix[0][1] = C_Yaw_S_Pitch * S_Roll - S_Yaw * C_Roll;
  Matrix[0][2] = C_Yaw_S_Pitch * C_Roll + S_Yaw * S_Roll;
  Matrix[0][3] = P1;
  Matrix[1][0] = S_Yaw * C_Pitch;
  Matrix[1][1] = S_Yaw_S_Pitch * S_Roll + C_Yaw * C_Roll;
  Matrix[1][2] = S_Yaw_S_Pitch * C_Roll - C_Yaw * S_Roll;
  Matrix[1][3] = P2;
  Matrix[2][0] = -S_Pitch;
  Matrix[2][1] = C_Pitch * S_Roll;
  Matrix[2][2] = C_Pitch * C_Roll;
  Matrix[2][3] = P3;
  Matrix[3][0] = Matrix[3][1] = Matrix[3][2] = 0;
  Matrix[3][3] = 1;
} /* ComputeT */


/**************************************************************************/
/* This routine computes the position factor between TI and TF            */
/**************************************************************************/
double ComputePosFact (DenType TI, DenType TF)
{
  register int i, j;
  double    Sum, Arm;

  for (i = 0, Sum = 0, Arm = 0; i < DenavitSize - 1; i++)
    for (j = 0; j < DenavitSize; j++)
      {
        Sum += pow ((TI[i][j] - TF[i][j]), 2);
        Arm += pow (TI[i][j], 2);
      } /* for */
  Sum = sqrt (Sum / Arm);
  return (Sum);
} /* ComputePosFact */


/**************************************************************************/
/* This routine fills the Y_P_R array with the Yaw_Pitch_Roll parameters  */
/* of the given frame TI. This parameters will be shown on the screen for */
/* different frames.                                                      */
/**************************************************************************/
void ComputeOrientation (DenType TI)
{
  Y_P_R [2] = atan2 (TI[2][1], TI[2][2]);  /* Roll  */
  Y_P_R [1] = atan2 (-TI[2][0], sqrt (pow (TI[0][0], 2) + pow (TI[1][0], 2)));
                                           /* Pitch */
  Y_P_R [0] = atan2 (TI[1][0], TI[0][0]);  /* Yaw   */

} /* ComputeOrientation */


/**************************************************************************/
/* This routine prints a message on the specified position of the screen  */
/* and if demanded clears it after 100ms.                                 */
/**************************************************************************/
void Message (char *S, int Col, int Clear)
{
  outtextxy (10, Col, S);
  delay (100);
  if (Clear)
    bar (10, Col, 250, Col+10);
} /* Message */


/**************************************************************************/
/* This routine prints L,M,N,P and Yaw_Pitch_Roll parameters of the given */
/* frame T at the given position. It also prints current joint angles in  */
/* its place.                                                             */
/**************************************************************************/
void PrintResult (DenType T, int Column)
{
  register int i, j;

  bar (5, Column, 309, Column+70);
  for (i = 0; i < DenavitSize; i++)
    {
      sprintf (Str, "%c ", Vector[i]);
      outtextxy (5, Column+i*10, Str);
      for (j = 0; j < DenavitSize - 1; j++)
        {
          sprintf (Str, "% 10.6lf", T[j][i]);
          outtextxy (20+j*100, Column+(i*10), Str);
        } /* for */
    } /* for */

  ComputeOrientation (T);
  for (i = 0; i < 3; i++)
    {
      sprintf (Str, "%s" , YPR[i]);
      outtextxy (40+i*100, Column+5*10, Str);
      sprintf (Str, "% 10.6lf", Y_P_R[i] * 180 / M_PI);
      outtextxy (20+i*100, Column+6*10, Str);
    } /* for */

  bar (5, 350, 620, 360);
  for (i = 0; i < MAXJOINT; i++)
    {
      sprintf (Str, "JOINT%d", i+1);
      outtextxy (40+i*100, 340, Str);
      sprintf (Str, "% 10.6lf", JointAngle[i] * 180 / M_PI);
      outtextxy (20+i*100, 350, Str);
    } /* for */

  bar (190, 380, 500, 390);
  PositionFactor = ComputePosFact (T_Initial, T_Final);
  sprintf (Str, "PositionFactor = % 10.6lf", PositionFactor);
  outtextxy (190, 380, Str);

  delay (20);
} /* PrintResult */


/***************************************************************************/
/* This routine checks if the distance between the two given frames TI and */
/* TF is less than ERRMAX; in other words if the two frames are converged. */
/***************************************************************************/
int Converge (DenType TI, DenType TF)
{
  register int i, j;
  double   Sum;

  for (i = 0, Sum = 0; i < DenavitSize - 1; i++)
    for (j = 0; j < DenavitSize; j++)
      Sum += pow ((TI[i][j] - TF[i][j]), 2);
  Sum = sqrt (Sum);
  return (Sum < ERRMAX);
} /* Converge */


/**************************************************************************/
/* After calling Initialization, until the convergence of the initial and */
/* final frames, if there is no other intermediate points for at most     */
/* MAXTRY times Interpolate is called to make the current and final frames*/
/* converge. In the presence of yet another intermediate point            */
/* Interpolate is called to make the current and next intermediate frams  */
/* converge, and then ResolveVibration is called to resolve any present   */
/* vibration or detect the absence of any more convergence.               */
/* If total convergence is gained or absence of convergence is detected   */
/* the main loop is exited and proper messages regarding convergence      */
/* the total number of steps, and the number of Vibration resolution steps*/
/* are printed in the message box.                                        */
/**************************************************************************/
void main (void)
{
  register int Convergence = 0, i;

  Initialization ();
  while (!Converge (T_Initial, T_Final))
    {
      if (!ComputeIntermediatePoint ())
        {
          for (i = 0; i < MAXTRY; i++)
            if (!(Convergence = Interpolate (T_Initial, T_Final)))
              break;
          break;
        } /* if */
      if (ResolveVibration (Convergence = Interpolate (T_Initial, T_Intermediate)))
        break;
    } /* while */
  if (Convergence)
    Message ("No Convergence", 420, 0);
  sprintf (Str, "Totally Done in %d Steps", Count);
  Message (Str, 430, 0);
  sprintf (Str, "With %d Vibration Resolution Steps", VibrateCount);
  Message (Str, 440, 0);
  getch ();
  closegraph();
} /* main */
