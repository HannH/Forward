#ifndef CGCORDINATE_H
#define CGCORDINATE_H
#include "DataStruct.h"

//const double PI = 3.1415926535897932384626433832795;
//const double DEG2RAD = PI / 180;
//const double RAD2DEG = 180 / PI;
//const double DOUBLE_0 = 1E-30;
//const double ITER_0 = 1E-15; //迭代收敛差值

#define DOUBLE_0 1E-30
#define ITER_0   1E-15

class Cg
{//你需要输入a 和 e的平方...长半轴a  第一偏心率e
	//长半轴，和第一偏心率
	//定义一个参考椭球需要两个参数a，e

	//WGS-84 a = 6378137, e2 = 0.0066943799013、
public:
	Cg(double _a, double _e2): a(_a), e2(_e2){}
public:
	double get_a() { return a; }
	double get_e2() { return e2; }

	void r2g(const double X, const double Y, const double Z, double &B, double &L, double &H);
	void g2r(const double B, const double L, const double H, double &X, double &Y, double &Z);
protected:
	double a;
	double e2;
};

void StaXYZ2BLH(Station3D st_XYZ, Station3DBLH &st_BLH, double a, double e2);
void StaBLH2XYZ(Station3DBLH st_BLH, Station3D &st_XYZ, double a, double e2);

void XYZ2BLH(SpacePoint3D sp_XYZ, SpacePoint3D_BLH &sp_BLH, double a, double e2);
void BLH2XYZ(SpacePoint3D_BLH sp_BLH, SpacePoint3D &sp_XYZ, double a, double e2);
#endif