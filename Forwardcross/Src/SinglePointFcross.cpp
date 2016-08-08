#include "SinglePointFcoss.h"
#include "math.h"

#include "Matrix.h"

//线性求解前方交会
CsinglePointForwardCross::CsinglePointForwardCross()
{
	memset(m_Qxx, 0, 9 * sizeof(double));
}

CsinglePointForwardCross::~CsinglePointForwardCross()
{
	if (!vpt2d.empty()) vpt2d.~vector();
	if (!vSt3d.empty()) vSt3d.~vector();
}

void CsinglePointForwardCross::CalRmatrix(Station3D st, double *R)
{
	double phi = st.phi, omega = st.omega, kappa = st.kappa;
	switch (st.rotate){
	case pok:{//p-o-k系统
			   double cp = cos(phi), sp = sin(phi);
			   double co = cos(omega), so = sin(omega);
			   double ck = cos(kappa), sk = sin(kappa);

			   R[0] = cp*ck - so*sp*sk; R[1] = -cp*sk - so*sp*ck; R[2] = -co*sp;
			   R[3] = co*sk; R[4] = co*ck; R[5] = -so;
			   R[6] = sp*ck + so*cp*sk; R[7] = -sp*sk + so*cp*ck; R[8] = co*cp;
	}
		break;
	case opk:{//-o-p-k系统
			   double co = cos(omega), so = sin(omega);
			   double cp = cos(phi), sp = sin(phi);
			   double ck = cos(kappa), sk = sin(kappa);

			   R[0] = cp * ck, R[1] = -cp * sk, R[2] = sp;
			   R[3] = so * sp * ck + co * sk, R[4] = -so * sp * sk + co * ck, R[5] = -so * cp;
			   R[6] = -co * sp * ck + so * sk, R[7] = co * sp * sk + so * ck, R[8] = co * cp;
	}
		break;
	default:
		fprintf(stderr, "转角系统有误！ 0 - p-o-k; 1 - o - p - k");
		break;
	}
}

bool CsinglePointForwardCross::DataPrepareWarning()
{
	//printf("正在进行线性求解前方交会!\n");
	if (vpt2d.size() != vSt3d.size()) {
		printf("光线数目不对！\n");
		false;
	}
	for (int i = 0; i < vSt3d.size(); i++) {
		if (vpt2d[i].imgID != vSt3d[i].imgID) {
			printf("像点和影像不对应！imgid = %d staId = %d\n", vpt2d[i].imgID, vSt3d[i].imgID);
			false;
		}

	}

	//printf("数据准备完成！\n");

	/*for (int i = 0; i < vSt3d.size(); i++) {
		printf("alpha = %lf beta = %lf, Xs = %lf, Ys = %lf, Zs = %lf\n",
		vpt2d[i].alpha * RAD2DEG, vpt2d[i].beta * RAD2DEG, vSt3d[i].Xs, vSt3d[i].Ys, vSt3d[i].Zs);
		}*/

	return true;
}

void CsinglePointForwardCross::ForwardCross()
{
	if(!DataPrepareWarning()) return;
	
	int equation_num = 2 * GetBundlLineNum(); //一条光线两个方程

	double *A   = new double[equation_num * 3]; 	memset(  A, 0, sizeof(double) * equation_num * 3);
	double *AT  = new double[3 * equation_num]; 	memset( AT, 0, sizeof(double) * equation_num * 3);
	double *L   = new double[equation_num * 1];     memset(  L, 0, sizeof(double) * equation_num * 1);
	double *X   = new double[3 * 1];                memset(  X, 0, sizeof(double) * 3 * 1);
	double *ATA = new double[3 * 3];                memset(ATA, 0, sizeof(double) * 3 * 3);
	double *ATL = new double[3 * 1];                memset(ATA, 0, sizeof(double) * 3 * 1);
	
	double l1, l2, l3, l4, l5, l6, lx, ly;
	double R[9];
	double Xs, Ys, Zs;
	double alpha_new, beta_new;
	int n = GetBundlLineNum();

	for (int i = 0, m = 0, j = 0; i < n; i++)	{
		CalRmatrix(vSt3d[i], R);

		Xs = vSt3d[i].Xs; Ys = vSt3d[i].Ys; Zs = vSt3d[i].Zs;
		alpha_new = tan(vpt2d[i].alpha); beta_new = tan(vpt2d[i].beta);

		l1 = R[0] + R[2] * alpha_new; l2 = R[3] + R[5] * alpha_new; l3 = R[6] + R[8] * alpha_new;
		l4 = R[1] + R[2] * beta_new ; l5 = R[4] + R[5] * beta_new ; l6 = R[7] + R[8] * beta_new ;

		lx = R[2] * alpha_new * Xs + R[5] * alpha_new * Ys + R[8] * alpha_new * Zs + R[0] * Xs + R[3] * Ys + R[6] * Zs;
		ly = R[2] * beta_new  * Xs + R[5] * beta_new  * Ys + R[8] * beta_new  * Zs + R[1] * Xs + R[4] * Ys + R[7] * Zs;

		A[m * 3 + 0]       = l1; A[m * 3 + 1]       = l2; A[m * 3 + 2]       = l3;
		A[(m + 1) * 3 + 0] = l4; A[(m + 1) * 3 + 1] = l5; A[(m + 1) * 3 + 2] = l6;	m += 2;
		L[j * 1 + 0]       = lx; L[(j + 1) * 1 + 0] = ly;                       	j += 2;//下一对同名点
	}

	CMatrix MatrixApp;
	MatrixApp.MatrixTranspose(A, AT, equation_num, 3);
	MatrixApp.MatrixMulti(AT, A, ATA, 3, equation_num, 3);
	MatrixApp.MatrixMulti(AT, L, ATL, 3, equation_num, 1);

	//printf("条件数：%lf\n", MatrixApp.Matrix_Condition(ATA, 3, 3)); 会改变矩阵的状态
	MatrixApp.MatrixInversion(ATA, 3);//ATA为协方差矩阵Qxx
	MatrixApp.MatrixMulti(ATA, ATL, X, 3, 3, 1);
	m_X = X[0]; m_Y = X[1]; m_Z = X[2];

	//double *V = new double[equation_num * 1];	memset(V, 0, sizeof(double)* equation_num);
	double *AX = new double[equation_num * 1];   memset(AX, 0, sizeof(double) * equation_num);
	MatrixApp.MatrixMulti(A, X, AX, equation_num, 3, 1);
	MatrixApp.FillMatrix_MINUS(AX, L, equation_num, 1, equation_num, 1, 0, 0); //V = AX - L ,结果放在V = AX

	double m0 = 0.0;//单位权方差
	for (int i = 0; i < equation_num; i++) {
		m0 += AX[i] * AX[i];
	}

	m0 = m0 / (equation_num - 3);

	//精度
	m_segma2 = m0;
	memcpy(m_Qxx, ATA, 9 * sizeof(double));

	/*m_DXX = m0 * ATA[0];
	m_DYY = m0 * ATA[4];
	m_DZZ = m0 * ATA[8];*/
	
	if(A)   delete[] A;
	if(AT)  delete[] AT;
	if (L)   delete[] L;
	if(X)   delete[] X;
	if(ATA) delete[] ATA;
	if(ATL) delete[] ATL;
	if(AX)  delete[] AX;
}


//迭代法求解前方交会
CintersinglePointForwardCross::CintersinglePointForwardCross()
{
	memset(m_Qxx, 0, 9 * sizeof(double));
}

CintersinglePointForwardCross::~CintersinglePointForwardCross()
{
	if (!m_vpt2d.empty()) m_vpt2d.~vector();
	if (!m_vst3d.empty()) m_vst3d.~vector();
}

void CintersinglePointForwardCross::CalRmatrix(Station3D st, double *R)
{
	double phi = st.phi, omega = st.omega, kappa = st.kappa;
	switch (st.rotate){
	case pok:{//p-o-k系统
				 double cp = cos(phi), sp = sin(phi);
				 double co = cos(omega), so = sin(omega);
				 double ck = cos(kappa), sk = sin(kappa);

				 R[0] = cp*ck - so*sp*sk; R[1] = -cp*sk - so*sp*ck; R[2] = -co*sp;
				 R[3] = co*sk; R[4] = co*ck; R[5] = -so;
				 R[6] = sp*ck + so*cp*sk; R[7] = -sp*sk + so*cp*ck; R[8] = co*cp;
	}
		break;
	case opk:{//-o-p-k系统
				 double co = cos(omega), so = sin(omega);
				 double cp = cos(phi), sp = sin(phi);
				 double ck = cos(kappa), sk = sin(kappa);

				 R[0] = cp * ck, R[1] = -cp * sk, R[2] = sp;
				 R[3] = so * sp * ck + co * sk, R[4] = -so * sp * sk + co * ck, R[5] = -so * cp;
				 R[6] = -co * sp * ck + so * sk, R[7] = co * sp * sk + so * ck, R[8] = co * cp;
	}
		break;
	default:
		fprintf(stderr, "转角系统有误！ 0 - p-o-k; 1 - o - p - k");
		break;
	}
}

bool CintersinglePointForwardCross::DataPrepareWarning()
{
	printf("正在进行迭代求解前方交会!\n");
	if (m_vpt2d.size() != m_vst3d.size()) {
		printf("光线数目不对！\n");
		false;
	}
	for (int i = 0; i < m_vst3d.size(); i++) {
		if (m_vpt2d[i].imgID != m_vst3d[i].imgID) {
			printf("像点和影像不对应！\n");
			false;
		}

	}

	printf("数据准备完成！\n");

	/*for (int i = 0; i < m_vst3d.size(); i++) {
		printf("alpha = %lf beta = %lf, Xs = %lf, Ys = %lf, Zs = %lf\n",
		m_vpt2d[i].alpha * RAD2DEG, m_vpt2d[i].beta * RAD2DEG, m_vst3d[i].Xs, m_vst3d[i].Ys, m_vst3d[i].Zs);
		}*/

	return true;
}

void CintersinglePointForwardCross::InterationForwardCross()
{
	if (!DataPrepareWarning()) return;

	int equation_num = 2 * GetBundlLineNum(); //一条光线两个方程

	double *A   = new double[equation_num * 3]; 	memset(A, 0, sizeof(double)* equation_num * 3);
	double *AT  = new double[3 * equation_num]; 	memset(AT, 0, sizeof(double)* equation_num * 3);
	double *L   = new double[equation_num * 1];     memset(L, 0, sizeof(double)* equation_num * 1);
	double *X   = new double[3 * 1];                memset(X, 0, sizeof(double)* 3 * 1);
	double *ATA = new double[3 * 3];                memset(ATA, 0, sizeof(double)* 3 * 3);
	double *ATL = new double[3 * 1];                memset(ATA, 0, sizeof(double)* 3 * 1);

	int ncount = 0; //迭代次数
	double R[9];
	double Xs, Ys, Zs;
	double X_new, Y_new, Z_new;
	double alpha_new, beta_new, alpha, beta;
	double l1, l2, l3, l4, l5, l6, lx, ly;

	int n = GetBundlLineNum();
	double X0 = m_X0, Y0 = m_Y0, Z0 = m_Z0; //初值
	CMatrix MatrixApp;
	do 
	{
		if (ncount > m_maxiter) break; //迭代次数限制
		ncount++;
		for (int i = 0, m = 0, j = 0; i < n; i++)	{
			CalRmatrix(m_vst3d[i], R);

			Xs = m_vst3d[i].Xs; Ys = m_vst3d[i].Ys; Zs = m_vst3d[i].Zs;
			
			X_new = R[0] * (X0 - Xs) + R[3] * (Y0 - Ys) + R[6] * (Z0 - Zs);
			Y_new = R[1] * (X0 - Xs) + R[4] * (Y0 - Ys) + R[7] * (Z0 - Zs);
			Z_new = R[2] * (X0 - Xs) + R[5] * (Y0 - Ys) + R[8] * (Z0 - Zs);

			alpha = atan(-X_new / Z_new); beta = atan(-Y_new / Z_new); //近似值
			alpha_new = tan(alpha); beta_new = tan(beta); //观测值正切

			l1 = -Z_new * (R[0] + R[2] * alpha_new) / (X_new * X_new + Z_new * Z_new);//dx/dX
			l2 = -Z_new * (R[3] + R[5] * alpha_new) / (X_new * X_new + Z_new * Z_new);//dx/dY
			l3 = -Z_new * (R[6] + R[8] * alpha_new) / (X_new * X_new + Z_new * Z_new);//dx/dZ
			l4 = -Z_new * (R[1] + R[2] * beta_new)  / (Y_new * Y_new + Z_new * Z_new);//dy/dX
			l5 = -Z_new * (R[4] + R[5] * beta_new)  / (Y_new * Y_new + Z_new * Z_new);//dy/dY
			l6 = -Z_new * (R[7] + R[8] * beta_new)  / (Y_new * Y_new + Z_new * Z_new);//dz/dZ
		
			lx = m_vpt2d[i].alpha - alpha; ly = m_vpt2d[i].beta - beta; //

			A[m * 3 + 0]       = l1; A[m * 3 + 1]       = l2; A[m * 3 + 2]       = l3;
			A[(m + 1) * 3 + 0] = l4; A[(m + 1) * 3 + 1] = l5; A[(m + 1) * 3 + 2] = l6;	m += 2;
			L[j * 1 + 0]       = lx; L[(j + 1) * 1 + 0] = ly;                       	j += 2;//下一对同名点
		}

		MatrixApp.MatrixTranspose(A, AT, equation_num, 3);
		MatrixApp.MatrixMulti(AT, A, ATA, 3, equation_num, 3);
		MatrixApp.MatrixMulti(AT, L, ATL, 3, equation_num, 1);

		//printf("条件数：%lf\n", MatrixApp.Matrix_Condition(ATA, 3, 3)); 会改变矩阵的状态
		MatrixApp.MatrixInversion(ATA, 3);//ATA为协方差矩阵Qxx
		MatrixApp.MatrixMulti(ATA, ATL, X, 3, 3, 1); //X为位置改正数

		X0 += X[0]; Y0 += X[1]; Z0 += X[2]; //改正， 迭代
	} while (X[0] > m_iter0 || X[1] > m_iter0 || X[2] > m_iter0 || ncount < m_miniter);

	printf("迭代次数： %d\n", ncount);
	m_X = X0; m_Y = Y0; m_Z = Z0;

	double *AX = new double[equation_num * 1];   memset(AX, 0, sizeof(double)* equation_num);
	MatrixApp.MatrixMulti(A, X, AX, equation_num, 3, 1);
	MatrixApp.FillMatrix_MINUS(AX, L, equation_num, 1, equation_num, 1, 0, 0); //V = AX - L ,结果放在V = AX

	printf("X:\n");
	for (int i = 0; i < 3; i++) printf("%lf ", X[i]);

	double m0 = 0.0;//单位权方差
	for (int i = 0; i < equation_num; i++) {
		m0 += AX[i] * AX[i];
	}

	m0 = m0 / (equation_num - 3);

	//精度
	m_m0 = m0;
	memcpy(m_Qxx, ATA, 9 * sizeof(double));

	if (A)   delete[] A;
	if (AT)  delete[] AT;
	if (L)   delete[] L;
	if (X)   delete[] X;
	if (ATA) delete[] ATA;
	if (ATL) delete[] ATL;
	if (AX)  delete[] AX;
}

