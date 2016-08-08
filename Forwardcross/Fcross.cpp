
#include "Src\SGCordinate.h"
#include "Src\Project.h"
#include "Src\Matrix.h"
#include "ExtraDef.h"
#include "Fcross.h"


CBundleAdjustment::CBundleAdjustment()
{

}
CBundleAdjustment:: ~CBundleAdjustment()
{

}

//该函数计算模型有待检查
void CBundleAdjustment::ForwardCross_11()
{
	int equation_num = 2 * m_bundlenum; //一条光线两个方程

	double *A = new double[equation_num * 3]; 	memset(A, 0, sizeof(double)* equation_num * 3);
	double *AT = new double[3 * equation_num]; 	memset(AT, 0, sizeof(double)* equation_num * 3);
	double *L = new double[equation_num * 1];     memset(L, 0, sizeof(double)* equation_num * 1);
	double *X = new double[3 * 1];                memset(X, 0, sizeof(double)* 3 * 1);
	double *ATA = new double[3 * 3];                memset(ATA, 0, sizeof(double)* 3 * 3);
	double *ATL = new double[3 * 1];                memset(ATA, 0, sizeof(double)* 3 * 1);

	double l1, l2, l3, l4, l5, l6, lx, ly;
	double R[9], RT[9], T0[9], K[9], T[9], C[3], S[3];
	double Xs, Ys, Zs;
	double alpha_new, beta_new;
	double alpha_new0, beta_new0;
	int n = m_bundlenum;

	CMatrix MatrixApp;
	for (int i = 0, m = 0, j = 0; i < n; i++)
	{
		//计算R， K, T
		switch (m_cmr[i].rotate) {
		case pok: POK2R(m_cmr[i].p, m_cmr[i].o, m_cmr[i].k, R); break;
		case opk: OPK2R(m_cmr[i].o, m_cmr[i].k, m_cmr[i].k, R); break;
		default: break;
		}

		//计算地心地固到轨道坐标系的旋转矩阵T0 = B * A
		//B*A*T + S - C才是想空间辅助坐标
		CalCoordinateTransfer(m_stBLH[i].L * DEG2RAD, m_stBLH[i].B * DEG2RAD, 0, T0);

		//计算旋转参数K = RT * T0
		MatrixApp.MatrixTranspose(R, RT, 3, 3);
		MatrixApp.MatrixMulti(RT, T0, K, 3, 3, 3);

		//计算平移参数T = RT * (S-C);
		C[0] = m_cmr[i].x0; C[1] = m_cmr[i].y0; C[2] = m_cmr[i].z0; //传感器到卫星的距离矢量
		Cg cg(WGS84_a, WGS84_e2);
		cg.g2r(m_stBLH[i].B, m_stBLH[i].L, m_stBLH[i].H, Xs, Ys, Zs);
		double dis = Xs * Xs + Ys * Ys + Zs * Zs;
		dis = sqrt(dis); //到地心距离
		S[0] = S[1] = 0;  S[2] = -dis; //到地心的位移矢量

		MatrixApp.FillMatrix_MINUS(S, C, 3, 1, 3, 1, 0, 0);
		MatrixApp.MatrixMulti(RT, S, T, 3, 3, 1);

		//[x, y, -f] = K * [X, Y, Z] + T
		alpha_new = tan(m_vpt2d[i].alpha); beta_new = tan(m_vpt2d[i].beta);
		alpha_new0 = tan(m_cmr[i].u);      beta_new0 = tan(m_cmr[i].v);

		l1 = K[0] + K[6] * (alpha_new - alpha_new0); 
		l2 = K[1] + K[7] * (alpha_new - alpha_new0); 
		l3 = K[2] + K[8] * (alpha_new - alpha_new0);
		l4 = K[3] + K[6] * (beta_new - beta_new0); 
		l5 = K[4] + K[7] * (beta_new - beta_new0); 
		l6 = K[5] + K[8] * (beta_new - beta_new0);

		lx = -T[0] - (alpha_new - alpha_new0) * T[2];
		ly = -T[1] - (beta_new - beta_new0);

		A[m * 3 + 0] = l1; A[m * 3 + 1] = l2; A[m * 3 + 2] = l3;
		A[(m + 1) * 3 + 0] = l4; A[(m + 1) * 3 + 1] = l5; A[(m + 1) * 3 + 2] = l6;	m += 2;
		L[j * 1 + 0] = lx; L[(j + 1) * 1 + 0] = ly;                       	j += 2;//下一对同名点
	}

	MatrixApp.MatrixTranspose(A, AT, equation_num, 3);
	MatrixApp.MatrixMulti(AT, A, ATA, 3, equation_num, 3);
	MatrixApp.MatrixMulti(AT, L, ATL, 3, equation_num, 1);

	//printf("条件数：%lf\n", MatrixApp.Matrix_Condition(ATA, 3, 3)); 会改变矩阵的状态
	MatrixApp.MatrixInversion(ATA, 3);//ATA为协方差矩阵Qxx
	MatrixApp.MatrixMulti(ATA, ATL, X, 3, 3, 1);

	m_vpt3d.X = X[0]; m_vpt3d.Y = X[1]; m_vpt3d.Z = X[2];

	double *AX = new double[equation_num * 1];   memset(AX, 0, sizeof(double)* equation_num);
	MatrixApp.MatrixMulti(A, X, AX, equation_num, 3, 1);
	MatrixApp.FillMatrix_MINUS(AX, L, equation_num, 1, equation_num, 1, 0, 0); //V = AX - L ,结果放在V = AX

	double m0 = 0.0;//单位权方差
	for (int i = 0; i < equation_num; i++) {
		m0 += AX[i] * AX[i];
	}

	m0 = m0 / (equation_num - 3);
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

void CBundleAdjustment::CalRmatrix(Station3D st, double *R)
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

void CBundleAdjustment::CalRmatrix_BLH(Station3DBLH st_BLH, CAMERA cmr, double *R)
{
	st_BLH.phi = cmr.p;
	st_BLH.omega = cmr.o;
	st_BLH.kappa = cmr.k;
	st_BLH.rotate = cmr.rotate;

	Station3D st;
	StaBLH2XYZ(st_BLH, st, WGS84_a, WGS84_e2);
	
	CalRmatrix(st, R);
}
#ifdef _DEBUG
#include <iostream>
#endif // _DEBUG

void CBundleAdjustment::ForwardCross()
{
	int equation_num = 2 * m_bundlenum; //一条光线两个方程

	double *A = new double[equation_num * 3]; 	memset(A, 0, sizeof(double)* equation_num * 3);
	double *AT = new double[3 * equation_num]; 	memset(AT, 0, sizeof(double)* equation_num * 3);
	double *L = new double[equation_num * 1];     memset(L, 0, sizeof(double)* equation_num * 1);
	double *X = new double[3 * 1];                memset(X, 0, sizeof(double)* 3 * 1);
	double *ATA = new double[3 * 3];                memset(ATA, 0, sizeof(double)* 3 * 3);
	double *ATL = new double[3 * 1];                memset(ATA, 0, sizeof(double)* 3 * 1);

	double l1, l2, l3, l4, l5, l6, lx, ly;
	double R[9];
	double Xs, Ys, Zs,x,y;
	double alpha_new, beta_new;
	int n = m_bundlenum;

	for (int i = 0, m = 0, j = 0; i < n; i++)	{
		CalRmatrix(m_cmr1[i],R);
		Xs = m_cmr1[i].Xs; Ys = m_cmr1[i].Ys; Zs = m_cmr1[i].Zs;
		x = m_vpt2d[i].x; y = m_vpt2d[i].y;
		double f = m_innerEle.Z, x0 = m_innerEle.X, y0 = m_innerEle.Y;
		l1 = f*R[0] + R[2] * (x - x0); l2 = f*R[3] + R[5] * (x - x0); l3 = f*R[6] + R[8] * (x - x0);
		l4 = f*R[1] + R[2] * (y - y0); l5 = f*R[4] + R[5] * (y - y0); l6 = f*R[7] + R[8] * (y - y0);

		lx = R[2] * (x - x0) * Xs + R[5] * (x - x0) * Ys + R[8] * (x - x0) * Zs + f*R[0] * Xs + f*R[3] * Ys + f*R[6] * Zs;
		ly = R[2] * (y - y0) * Xs + R[5] * (y - y0) * Ys + R[8] * (y - y0) * Zs + f*R[1] * Xs + f*R[4] * Ys + f*R[7] * Zs;

		A[m * 3 + 0] = l1; A[m * 3 + 1] = l2; A[m * 3 + 2] = l3;
		A[(m + 1) * 3 + 0] = l4; A[(m + 1) * 3 + 1] = l5; A[(m + 1) * 3 + 2] = l6;	m += 2;
		L[j] = lx; L[j + 1] = ly;                       	j += 2;//下一对同名点
	}

	CMatrix MatrixApp;
	MatrixApp.MatrixTranspose(A, AT, equation_num, 3);
	MatrixApp.MatrixMulti(AT, A, ATA, 3, equation_num, 3);
	MatrixApp.MatrixMulti(AT, L, ATL, 3, equation_num, 1);

	//printf("条件数：%lf\n", MatrixApp.Matrix_Condition(ATA, 3, 3)); 会改变矩阵的状态
	MatrixApp.MatrixInversion(ATA, 3);//ATA为协方差矩阵Qxx
	MatrixApp.MatrixMulti(ATA, ATL, X, 3, 3, 1);
	m_vpt3d.X = X[0]; m_vpt3d.Y = X[1]; m_vpt3d.Z = X[2];

	//double *V = new double[equation_num * 1];	memset(V, 0, sizeof(double)* equation_num);
	double *AX = new double[equation_num * 1];   memset(AX, 0, sizeof(double)* equation_num);
	MatrixApp.MatrixMulti(A, X, AX, equation_num, 3, 1);
	MatrixApp.FillMatrix_MINUS(AX, L, equation_num, 1, equation_num, 1, 0, 0); //V = AX - L ,结果放在V = AX

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