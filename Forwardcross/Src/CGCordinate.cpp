#include "SGCordinate.h"

#include "math.h"
#include "Project.h"

void Cg::g2r(const double B, const double L, const double H, double &X, double &Y, double &Z)
{
	double radB = B*DEG2RAD;
	double radL = L*DEG2RAD;

	double sinB = sin(radB);
	double cosB = cos(radB);

	double N = a / sqrt(1 - e2*sinB*sinB);
	double NHcosB = (N + H) * cosB;

	X = NHcosB * cos(radL);
	Y = NHcosB * sin(radL);
	Z = (N*(1 - e2) + H) * sinB;
}


void Cg::r2g(const double X, const double Y, const double Z, double &B, double &L, double &H)
{
	double r = sqrt(X*X + Y*Y);

	if (abs(r) < DOUBLE_0) //极点处的特殊情况
	{
		L = 0; //L可以设置为任意值
		if (Z > 0)
		{
			B = 90;
			H = Z - a*sqrt(1 - e2); //大地测量学6-36的简化
		}
		else
		{
			B = -90;
			H = -Z - a*sqrt(1 - e2);
		}
	}
	else
	{
		//1.L
		if (abs(X) < DOUBLE_0) //L=90或-90时
		{
			L = asin(Y / r) * RAD2DEG;  // 可能存在符号问题
		}
		else
		{
			L = atan2(Y, X) * RAD2DEG;
		}
		//////////////////////////////////////////////////////////////////////////
		double t0 = Z / r;
		double c = a / sqrt(1 - e2);
		double P = c*e2 / r;
		double k = 1 / (1 - e2);

		double tanB0 = t0;
		double tanB1 = t0;

		do
		{
			tanB0 = tanB1;
			tanB1 = t0 + P*tanB0 / sqrt(k + tanB0*tanB0);
		} while (abs(tanB1 - tanB0) > ITER_0);

		double radB = atan(tanB1);
		B = radB * RAD2DEG;
		//3.H
		double sinB = sin(radB);
		double N = a / sqrt(1 - e2*sinB*sinB);
		H = r / cos(radB) - N;
	}
}

//将卫星大地坐标系下的摄站位置转化为经纬高
void StaXYZ2BLH(Station3D st_XYZ, Station3DBLH &st_BLH, double a, double e2)
{
	double B, L, H;
	Cg cg(a, e2);
	cg.r2g(st_XYZ.Xs, st_XYZ.Ys, st_XYZ.Zs, B, L, H);
	st_BLH.imgID = st_XYZ.imgID;
	st_BLH.B = B; st_BLH.L = L; st_BLH.H = H;
	st_BLH.omega = st_XYZ.omega; st_BLH.phi = st_XYZ.phi; st_BLH.kappa = st_XYZ.kappa;
	st_BLH.rotate = st_XYZ.rotate;
}
//将卫星经纬高的摄站位置转化为大地坐标系的坐标
void StaBLH2XYZ(Station3DBLH st_BLH, Station3D &st_XYZ, double a, double e2)
{
	double X, Y, Z;
	Cg cg(a, e2);
	cg.g2r(st_BLH.B, st_BLH.L, st_BLH.H, X, Y, Z);
	st_XYZ.imgID = st_BLH.imgID;
	st_XYZ.Xs = X; st_XYZ.Ys = Y; st_XYZ.Zs = Z;
	st_XYZ.omega = st_BLH.omega; st_XYZ.phi = st_BLH.phi; st_XYZ.kappa = st_BLH.kappa;
	st_XYZ.rotate = st_BLH.rotate;

	double p, o, k;

	double R[9];
	switch (st_XYZ.rotate){
	case pok: {
				  POK2R(st_XYZ.phi, st_XYZ.omega, st_XYZ.kappa, R);
				  CalPOK(st_BLH.L * DEG2RAD, st_BLH.B * DEG2RAD, 0, R, p, o, k);
				  st_XYZ.omega = o; st_XYZ.phi = p; st_XYZ.kappa = k;
	}break;
	case opk: {
				  OPK2R(st_XYZ.omega, st_XYZ.phi, st_XYZ.kappa, R);
				  CalOPK( st_BLH.L * DEG2RAD, st_BLH.B * DEG2RAD, 0,R, o, p, k);
				  st_XYZ.omega = o; st_XYZ.phi = p; st_XYZ.kappa = k;
	}break;
	default: break;
	}
}

void XYZ2BLH(SpacePoint3D sp_XYZ, SpacePoint3D_BLH &sp_BLH, double a, double e2)
{
	Cg cg(a, e2);
	cg.r2g(sp_XYZ.X, sp_XYZ.Y, sp_XYZ.Z, sp_BLH.B, sp_BLH.L, sp_BLH.H);
}
void BLH2XYZ(SpacePoint3D_BLH sp_BLH, SpacePoint3D &sp_XYZ, double a, double e2)
{
	Cg cg(a, e2);
	cg.g2r(sp_BLH.B, sp_BLH.L, sp_BLH.H, sp_XYZ.X, sp_XYZ.Y, sp_XYZ.Z);
}