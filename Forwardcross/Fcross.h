#ifndef FCROSS_H
#define FCROSS_H
//前方交会

#include "Src/DataStruct.h"

#pragma warning (disable : 4018)

class CBundleAdjustment
{
public:
	CBundleAdjustment();
	~CBundleAdjustment();
	
	void GetSpacePoint3D(SpacePoint3D &vpt3d) {	vpt3d = m_vpt3d; }
	void SetParameters(vector<Point2D> vpt2d, vector<CAMERA> cmr, vector<Station3DBLH> stBLH, int n) {
		SetPoint2D(vpt2d); SetCmrs(cmr); SetSatelliteLocation(stBLH); SetBundleNum(n);
	}
	void SetParameters(vector<Point2D> vpt2d, vector<Station3D> cmr, vector<SpacePoint3D> stPoint, int n) {
		SetPoint2D(vpt2d); SetCmrs(cmr); SetSatelliteLocation(stPoint); SetBundleNum(n);
	}
	void RunForwardCross() { ForwardCross(); }
	void SetInnerEle(int x0, int y0, int f) { m_innerEle.X = x0, m_innerEle.Y = y0, m_innerEle.Z = f; }

private:
	void SetPoint2D(vector<Point2D> vpt2d) { m_vpt2d.assign(vpt2d.begin(), vpt2d.end()); }
	void SetPoint2D(Point2D *vpt2d, int n) {
		for (int i = 0; i < n; i++) {
			m_vpt2d.push_back(vpt2d[i]);
		}
	}
	void SetCmrs(vector<CAMERA> cmr) { for (int i = 0; i < cmr.size(); i++) m_cmr.push_back(cmr[i]); }
	void SetCmrs(vector<Station3D> cmr) { for (int i = 0; i < cmr.size(); i++) m_cmr1.push_back(cmr[i]); }
	void SetSatelliteLocation(vector<Station3DBLH> stBLH) { for (int i = 0; i < stBLH.size(); i++) m_stBLH.push_back(stBLH[i]); }
	void SetSatelliteLocation(vector<SpacePoint3D> st) { for (int i = 0; i < st.size(); i++) m_st1.push_back(st[i]); }
	void SetBundleNum(int n) { m_bundlenum = n; }
	void ForwardCross_11();   //线性前方交会

	void ForwardCross();
	void CalRmatrix(Station3D st, double *R);//旋转矩阵
	void CalRmatrix_BLH(Station3DBLH st_BLH, CAMERA cmr, double *R);

	vector<Point2D> m_vpt2d;
	SpacePoint3D    m_vpt3d;
	vector<CAMERA>  m_cmr;
	//-----添加处
	vector<Station3D> m_cmr1;
	vector<SpacePoint3D> m_st1;
	SpacePoint3D m_innerEle;
	//--------完
	vector<Station3DBLH> m_stBLH;
	int m_bundlenum;  //光线数目
	double m_m0;
	double m_Qxx[9];  //协方差矩阵
};
#endif