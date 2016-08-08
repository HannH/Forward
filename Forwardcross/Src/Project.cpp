#include "Project.h"

#include "math.h"

#include "Matrix.h"

double roundx(double x) {
	double result;
	if (x > 1) result = 1;
	else if (x < -1) result = -1;
	else{ result = x; }
	return result;
}

Point2D invproject(Station3D st, SpacePoint3D pt3d,bool isPlane)
{
	double R[9];
	Point2D pt;

	CalRotate(st, R);

	CMatrix matrixapp;
	printf("R:\n");
	matrixapp.MatrixPrint(R, 3, 3);

	double Xs = st.Xs , Ys = st.Ys , Zs = st.Zs ;
	double X  = pt3d.X, Y  = pt3d.Y, Z  = pt3d.Z;

	double X_new = R[0] * (X - Xs) + R[3] * (Y - Ys) + R[6] * (Z - Zs);
	double Y_new = R[1] * (X - Xs) + R[4] * (Y - Ys) + R[7] * (Z - Zs);
	double Z_new = R[2] * (X - Xs) + R[5] * (Y - Ys) + R[8] * (Z - Zs);

	double x, y;
	if (abs(Z_new) <= 1e-15) {
		x = -atan2(X_new, Z_new);
		y = -atan2(Y_new, Z_new);
		return pt;
	}
	double f = -0.1;
	x = -f*X_new / Z_new;
	y = -f*Y_new / Z_new;
	if (isPlane == true) return Point2D(x, y);

	pt.alpha = atan(x); //东西视场角只在：-PI/2~PI/2， atan2...
	pt.beta  = atan(y);
	pt.imgID = st.imgID;

	return pt;
}

void CalRotate(Station3D st, double *R)
{
	double phi = st.phi, omega = st.omega, kappa = st.kappa;
	switch (st.rotate) {
	case pok:{//p-o-k系统
			   double cp = cos(phi)  , sp = sin(phi)  ;
			   double co = cos(omega), so = sin(omega);
			   double ck = cos(kappa), sk = sin(kappa);

			   R[0] = cp*ck - so*sp*sk; R[1] = -cp*sk - so*sp*ck; R[2] = -co*sp;
			   R[3] = co*sk           ; R[4] =  co*ck           ; R[5] = -so   ;
			   R[6] = sp*ck + so*cp*sk; R[7] = -sp*sk + so*cp*ck; R[8] =  co*cp;
	}
		break;
	case opk:{//-o-p-k系统
			   double co = cos(omega), so = sin(omega);
			   double cp = cos(phi)  , sp = sin(phi)  ;
			   double ck = cos(kappa), sk = sin(kappa);

			   R[0] =  cp * ck               , R[1] = -cp * sk               , R[2] =  sp     ;
			   R[3] =  so * sp * ck + co * sk, R[4] = -so * sp * sk + co * ck, R[5] = -so * cp;
			   R[6] = -co * sp * ck + so * sk, R[7] =  co * sp * sk + so * ck, R[8] =  co * cp;
	}
		break;
	default:
		fprintf(stderr, "转角系统有误！ 0 - p-o-k; 1 - o - p - k");
		break;
	}
}

void POK2R(const double p, const double o, const double k, double *R)
{
	double cp = cos(p), sp = sin(p);
	double co = cos(o), so = sin(o);
	double ck = cos(k), sk = sin(k);

	R[0] = cp*ck - so*sp*sk; R[1] = -cp*sk - so*sp*ck; R[2] = -co*sp;
	R[3] = co*sk; R[4] = co*ck; R[5] = -so;
	R[6] = sp*ck + so*cp*sk; R[7] = -sp*sk + so*cp*ck; R[8] = co*cp;
}

void OPK2R(const double o, const double p, const double k, double *R)
{
	double co = cos(o), so = sin(o);
	double cp = cos(p), sp = sin(p);
	double ck = cos(k), sk = sin(k);

	R[0] = cp * ck, R[1] = -cp * sk, R[2] = sp;
	R[3] = so * sp * ck + co * sk, R[4] = -so * sp * sk + co * ck, R[5] = -so * cp;
	R[6] = -co * sp * ck + so * sk, R[7] = co * sp * sk + so * ck, R[8] = co * cp;
}

//计算大地坐标系到轨道坐标系的旋转矩阵 T1 = A * T, T2 = B * T1, A = R(u)*R(i)*R(s)
//     s  - 升交点赤经   
//     i  - 轨道倾角     
//     u  - 纬度幅度（轨道近地点幅角+轨道真近点角）  逆时针为正   
//|      | 0  1  0 |            | cos(s) -sin(s)  0 |           |1     0         0|           | cos(u) -sin(u) 0 |
//|   B= | 0  0 -1 |  R(s)(Z) = | sin(s)  cos(s)  0 | R(i)(X) = |0  cos(i) -sin(i)| R(u)(Z) = | sin(u)  cos(u) 0 |
//|      |-1  0  0 |            |      0       0  1 |           |0  sin(i)  cos(i)|           |      0      0  1 |
//////////////////////////////////////////////////////////////////////////
void CalCoordinateTransfer(const double s, const double i, const double u, double *T)
{//2016-01-04,xxw
	double A[9];
	A[0] =  cos(u) * cos(s) - sin(u) * cos(i) * sin(s);
	A[1] = -cos(u) * sin(s) - sin(u) * cos(i) * cos(s);
	A[2] =  sin(u) * sin(i);
	A[3] =  sin(u) * cos(s) + cos(u) * cos(i) * sin(s);
	A[4] = -sin(u) * sin(s) + cos(u) * cos(i) * cos(s);
	A[5] = -cos(u) * sin(i);
	A[6] =  sin(i) * sin(s);
	A[7] =  sin(i) * cos(s);
	A[8] =  cos(i);

	//T = B*A
	T[0] =  A[3]; T[1] =  A[4]; T[2] =  A[5];
	T[3] = -A[6]; T[4] = -A[7]; T[5] = -A[8];
	T[6] = -A[0]; T[7] = -A[1]; T[8] = -A[2];
}

//////////////////////////////////////////////////////////////////////////
// 通过逆转化矩阵，计算旋转矩阵，并计算p-o-k的值
//     s  - 升交点赤经   (静止卫星为经度)
//     i  - 轨道倾角     
//     u  - 纬度幅度（轨道近地点幅角 + 轨道真近点角）
//     R0 - 卫星本体坐标系下的旋转矩阵
// R = inv(B*Ru*Ri*Rs) = inv(Rs) * inv(R(i)) * inv(Ru) * inv(B)
//         | 0 0 -1 |         |  cos(s)   sin(s) 0|           |1     0         0|             | cos(u)  sin(u)   0|
//inv(B) = | 1 0  0 | inv(Rs)=| -sin(s)   cos(s) 0| inv(R(i))=|0  cos(i)  sin(i)| inv(R(u)) = |-sin(u)  cos(u)   0|
//         | 0 -1 0 |         |     0       0    1|           |0 -sin(i)  cos(i)|             |     0        0   1|
void CalPOK(const double s, const double i, const double u, const double *R0, double &p, double &o, double &k)
{
	//double R[9];
	//R[0] =  cos(s) * cos(u) - sin(s) * cos(i) * sin(u);
	//R[1] = -cos(s) * sin(u) - sin(s) * cos(i) * cos(u);
	//R[2] =  sin(s) * sin(i);
	//R[3] =  sin(s) * cos(u) + cos(s) * cos(i) * sin(u);
	//R[4] = -sin(s) * sin(u) + cos(s) * cos(i) * cos(u);
	//R[5] = -cos(s) * sin(i);
	//R[6] =  sin(i) * sin(u);
	//R[7] =  sin(i) * cos(u);
	//R[8] =  cos(i);

	////R = R*inv(B)
	//double T[9]; //旋转矩阵
	//T[0] = R[1]; T[1] = -R[2]; T[2] = -R[0];
	//T[3] = R[4]; T[4] = -R[5]; T[5] = -R[3];
	//T[6] = R[7]; T[7] = -R[8]; T[8] = -R[6];

	double T[9];
	CalCoordinateTransfer(s, i, u, T);
	double Tinv[9];
	memcpy(Tinv, T, 9 * sizeof(double));

	CMatrix matrixapp;
	//matrixapp.MatrixPrint(T, 3, 3);
	matrixapp.MatrixInversion_General(Tinv, 3);

	double R[9];
	memcpy(R, R0, 9 * sizeof(double));

	matrixapp.MatrixMulti(Tinv, R, T, 3, 3, 3);

	//double R0_tmp[9];
	//memcpy(R0_tmp, R0, 9 * sizeof(double));

	//CMatrix MatrixApp;
	//MatrixApp.MatrixInversion(R0_tmp, 3);

	//double T_tmp[9];
	//MatrixApp.MatrixMulti(T, R0_tmp, T_tmp, 3, 3, 3); // 真实的旋转矩阵
	//memcpy(T, T_tmp, 9 * sizeof(double));

	o = -asin(roundx(T[5]));
	p = -asin(roundx(T[2] / cos(o)));
	k =  asin(roundx(T[3] / cos(o)));
	double k2 = acos(roundx(T[4] / cos(o)));

	p  = p  * RAD2DEG;
	o  = o  * RAD2DEG;
	k  = k  * RAD2DEG;
	k2 = k2 * RAD2DEG;

	// kappa range from -180~180
	if (k2 < 90){ ; }
	else if (k > 0){ k = k2; }
	else{ k = -k2; }

	//printf("p = %lf o = %lf k = %lf\n", p, o, k);

	p = p * DEG2RAD;//clockwise
	o = o * DEG2RAD;
	k = k * DEG2RAD;
}

void CalOPK(const double s, const double i, const double u, const double *R0, double &o, double &p, double &k)
{
	//double R[9];
	//R[0] =  cos(s) * cos(u) - sin(s) * cos(i) * sin(u);
	//R[1] = -cos(s) * sin(u) - sin(s) * cos(i) * cos(u);
	//R[2] =  sin(s) * sin(i);
	//R[3] =  sin(s) * cos(u) + cos(s) * cos(i) * sin(u);
	//R[4] = -sin(s) * sin(u) + cos(s) * cos(i) * cos(u);
	//R[5] = -cos(s) * sin(i);
	//R[6] =  sin(i) * sin(u);
	//R[7] =  sin(i) * cos(u);
	//R[8] =  cos(i);

	////R = R*inv(B)
	//double T[9]; //旋转矩阵
	//T[0] = R[1]; T[1] = -R[2]; T[2] = -R[0];
	//T[3] = R[4]; T[4] = -R[5]; T[5] = -R[3];
	//T[6] = R[7]; T[7] = -R[8]; T[8] = -R[6];

	//double R0_tmp[9];
	//memcpy(R0_tmp, R0, 9 * sizeof(double));

	//CMatrix MatrixApp;
	//MatrixApp.MatrixInversion(R0_tmp, 3);

	//double T_tmp[9];
	//MatrixApp.MatrixMulti(T, R0_tmp, T_tmp, 3, 3, 3); // 真实的旋转矩阵
	//memcpy(T, T_tmp, 9 * sizeof(double));

	double T[9];
	CalCoordinateTransfer(s, i, u, T);
	double Tinv[9];
	memcpy(Tinv, T, 9 * sizeof(double));

	CMatrix matrixapp;
	//matrixapp.MatrixPrint(T, 3, 3);
	matrixapp.MatrixInversion_General(Tinv, 3);

	double R[9];
	memcpy(R, R0, 9 * sizeof(double));

	matrixapp.MatrixMulti(Tinv, R, T, 3, 3, 3);

	p =  asin(roundx(T[2]));
	o = -asin(roundx(T[5] / cos(p)));
	k = -asin(roundx(T[1] / cos(p)));
	double k2 = acos(roundx(T[0] / cos(p)));

	p  = p  * RAD2DEG;
	o  = o  * RAD2DEG;
	k  = k  * RAD2DEG;
	k2 = k2 * RAD2DEG;

	if (k2 < 90){ ; }
	else if (k > 0) { k = k2; }
	else{ k = -k2; }

	p = p * DEG2RAD;//counterclockwise
	o = o * DEG2RAD;
	k = k * DEG2RAD;
}