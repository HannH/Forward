#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#include <stdlib.h>
#include <vector>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define DEG2RAD (PI/180)
#define RAD2DEG (180/PI)

#define WGS84_a 6378.137
#define WGS84_e2 0.0066943799013

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define  DOUBLE_MIN  99999999999999.0
#define  DOUBLE_MAX -99999999999999.0

//#define  pok 0
//#define  opk 1

enum RotateSym{
	pok,
	opk
};

enum RandomSym{
	Normal_Distribution,
	Uniform_Distribution
};

enum FcrossMethod{
	ForwardCross_Linear,
	ForwardCross_Inter
};
struct Point2D{//二维平面点
	int imgID;      //影像编号
	double alpha;   //东西方向,弧度,x
	double beta;    //南北方向,弧度,y
	double &x = alpha;
	double &y = beta;
	Point2D(double _a, double _b) : alpha(_a), beta(_b){}
	Point2D(int _ID, double _a, double _b) : imgID(_ID), alpha(_a), beta(_b){}

	Point2D & operator = (Point2D &pt) {
		imgID = pt.imgID; alpha = pt.alpha; beta = pt.beta;
		return *this;
	}
		Point2D(){
			imgID = -1;
			alpha = beta = 0.0;
		}
};

struct Station3D{//摄站点
	int imgID;  //影像编号
	double Xs, Ys, Zs;
	double omega, phi, kappa; 
	RotateSym rotate; //转角方式
	Station3D(double _Xs, double _Ys, double _Zs, double _o, double _p, double _k, RotateSym _r) :
		Xs(_Xs), Ys(_Ys), Zs(_Zs), omega(_o), phi(_p), kappa(_k), rotate(_r){
		//rotate = pok;
	}
	Station3D(int _ID, double _Xs, double _Ys, double _Zs, double _o, double _p, double _k, RotateSym _r) : imgID(_ID), 
		Xs(_Xs), Ys(_Ys), Zs(_Zs), omega(_o), phi(_p), kappa(_k), rotate(_r){
		//rotate = pok;
	}
	Station3D(){
		imgID = -1;
		Xs = Ys = Zs = 0.0;
		omega = phi = kappa = 0.0;
		rotate = pok;
	}
};

struct Station3DBLH{
	int imgID;  //影像编号
	double B, L, H;
	double omega, phi, kappa;//标称图的姿态
	RotateSym rotate; //转角方式
	Station3DBLH(double _B, double _L, double _H, double _o, double _p, double _k, RotateSym _r) :
		B(_B), L(_L), H(_H), omega(_o), phi(_p), kappa(_k), rotate(_r){
		//rotate = pok;
	}
	Station3DBLH(int _ID, double _B, double _L, double _H, double _o, double _p, double _k, RotateSym _r) : imgID(_ID),
		B(_B), L(_L), H(_H), omega(_o), phi(_p), kappa(_k), rotate(_r){
		//rotate = pok;
	}
	Station3DBLH(){
		imgID = -1;
		B = L = H = 0.0;
		omega = phi = kappa = 0.0;
		rotate = pok;
	}
};

struct SpacePoint3D{//目标空间点
	int spaceID;
	double X, Y, Z;
	SpacePoint3D(double _X, double _Y, double _Z) : X(_X), Y(_Y), Z(_Z){}
	SpacePoint3D(int _spaceID, double _X, double _Y, double _Z) : spaceID(_spaceID), X(_X), Y(_Y), Z(_Z){}
	SpacePoint3D & operator = (SpacePoint3D sp) {
		spaceID = sp.spaceID; X = sp.X; Y = sp.Y; Z = sp.Z; return *this;
	}
	SpacePoint3D(){
		spaceID = -1;
		X = Y = Z = 0.0;
	}
};

struct SpacePoint3D_BLH{
	int spaceID;
	double B, L, H;
	SpacePoint3D_BLH(double _B, double _L, double _H) : B(_B), L(_L), H(_H){}
	SpacePoint3D_BLH(int _spaceID, double _B, double _L, double _H) : spaceID(_spaceID), B(_B), L(_L), H(_H){}
	SpacePoint3D_BLH & operator = (SpacePoint3D_BLH sp) {
		spaceID = sp.spaceID; B = sp.B; L = sp.L; H = sp.H; return *this;
	}
	SpacePoint3D_BLH(){
		spaceID = -1;
		B = L = H = 0.0;
	}
};

struct BundleLine{//摄影光线
	Point2D pt;
	Station3D St3D;
};

//wang xiang
struct WxVirtualPoint{
	int ID;
	double time;
	SpacePoint3D pt3d; //大地坐标
	Point2D pl, pr; //左右像点
};

#endif