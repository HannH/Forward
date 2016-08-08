#ifndef SINGLEPOINTFCROSS_H
#define SINGLEPOINTFCROSS_H

#include "DataStruct.h"

//前方交会线性求解
class CsinglePointForwardCross{
public:
	CsinglePointForwardCross();
	~CsinglePointForwardCross();

	CsinglePointForwardCross(vector<Point2D> _vpt, vector<Station3D> _vst) : vpt2d(_vpt), vSt3d(_vst){}

	int GetBundlLineNum(){ return vSt3d.size(); } //光线数目
	void GetSpacePoint3D(SpacePoint3D &sp) { sp.X = m_X; sp.Y = m_Y; sp.Z = m_Z; }
	void GetSpacePoint3D_DD(double &Dxx, double &Dyy, double &Dzz) {
		Dxx = m_segma2 * m_Qxx[0]; Dyy = m_segma2 * m_Qxx[4]; Dzz = m_segma2 * m_Qxx[8]; }
	void GetSpacePoint3D_Qxx(double *Qxx) { memcpy(Qxx, m_Qxx, 9 * sizeof(double)); }
	void GetSpacePoint3D_segma2(double &segma2){ segma2 = m_segma2; }

	void RunSinglePointForwardCross(){ ForwardCross(); }

private:
	
	void CalRmatrix(Station3D st, double *R);//旋转矩阵
	void ForwardCross(); //前方交会
	bool DataPrepareWarning(); //准备光线

private:
	vector<Point2D> vpt2d; //二维平面点
	vector<Station3D> vSt3d; //摄站点
	double m_X;
	double m_Y;
	double m_Z;

	//精度
	double m_DXX;
	double m_DYY;
	double m_DZZ;
	double m_segma2; //单位权中误差

	double m_Qxx[9];  //协方差矩阵
};

//前方交会迭代求解
class CintersinglePointForwardCross{
public:
	CintersinglePointForwardCross();
	~CintersinglePointForwardCross();
	CintersinglePointForwardCross(vector<Point2D> _vpt, vector<Station3D> _vst, double _X0, double _Y0, double _Z0, double _iter0) :
		m_vpt2d(_vpt), m_vst3d(_vst), m_X0(_X0), m_Y0(_Y0), m_Z0(_Z0), m_iter0(_iter0)	{
		m_maxiter = 150; m_miniter = 3;
	}
public:
	int GetBundlLineNum(){ return m_vst3d.size(); } //光线数目
	void GetSpacePoint3D(SpacePoint3D &spt){ spt.X = m_X; spt.Y = m_Y; spt.Z = m_Z; }
	void GetSpacePoint3D_m0(double &m0) { m0 = m_m0; } //单位权方差
	void GetSpacePoint3D_Qxx(double *Qxx){ memcpy(Qxx, m_Qxx, 9 * sizeof(double)); }//协方差矩阵
	void GetSpacePoint3D_Dxx(double &Dxx, double &Dyy, double &Dzz){ //方差
		Dxx = m_m0 * m_Qxx[0]; Dyy = m_m0 * m_Qxx[4]; Dzz = m_m0 * m_Qxx[8]; }
	void RunintersinglePointForwardCross(){ InterationForwardCross(); }
private:
	void CalRmatrix(Station3D st, double *R);//旋转矩阵
	void InterationForwardCross(); //迭代法前方交会
	bool DataPrepareWarning(); //准备光线
private:
	//迭代初值
	double m_X0; //迭代求解的初值
	double m_Y0;
	double m_Z0;
	double m_iter0; //迭代精度（和XYZ一致）
	double m_maxiter; //最大迭代次数
	double m_miniter; //最小迭代次数

	//结果值
	double m_X;
	double m_Y;
	double m_Z;

	//单位权方差, 协方差矩阵
	double m_m0;
	double m_Qxx[9];

	vector<Point2D> m_vpt2d;
	vector<Station3D> m_vst3d;
};

#endif