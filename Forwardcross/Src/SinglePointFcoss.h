#ifndef SINGLEPOINTFCROSS_H
#define SINGLEPOINTFCROSS_H

#include "DataStruct.h"

//ǰ�������������
class CsinglePointForwardCross{
public:
	CsinglePointForwardCross();
	~CsinglePointForwardCross();

	CsinglePointForwardCross(vector<Point2D> _vpt, vector<Station3D> _vst) : vpt2d(_vpt), vSt3d(_vst){}

	int GetBundlLineNum(){ return vSt3d.size(); } //������Ŀ
	void GetSpacePoint3D(SpacePoint3D &sp) { sp.X = m_X; sp.Y = m_Y; sp.Z = m_Z; }
	void GetSpacePoint3D_DD(double &Dxx, double &Dyy, double &Dzz) {
		Dxx = m_segma2 * m_Qxx[0]; Dyy = m_segma2 * m_Qxx[4]; Dzz = m_segma2 * m_Qxx[8]; }
	void GetSpacePoint3D_Qxx(double *Qxx) { memcpy(Qxx, m_Qxx, 9 * sizeof(double)); }
	void GetSpacePoint3D_segma2(double &segma2){ segma2 = m_segma2; }

	void RunSinglePointForwardCross(){ ForwardCross(); }

private:
	
	void CalRmatrix(Station3D st, double *R);//��ת����
	void ForwardCross(); //ǰ������
	bool DataPrepareWarning(); //׼������

private:
	vector<Point2D> vpt2d; //��άƽ���
	vector<Station3D> vSt3d; //��վ��
	double m_X;
	double m_Y;
	double m_Z;

	//����
	double m_DXX;
	double m_DYY;
	double m_DZZ;
	double m_segma2; //��λȨ�����

	double m_Qxx[9];  //Э�������
};

//ǰ������������
class CintersinglePointForwardCross{
public:
	CintersinglePointForwardCross();
	~CintersinglePointForwardCross();
	CintersinglePointForwardCross(vector<Point2D> _vpt, vector<Station3D> _vst, double _X0, double _Y0, double _Z0, double _iter0) :
		m_vpt2d(_vpt), m_vst3d(_vst), m_X0(_X0), m_Y0(_Y0), m_Z0(_Z0), m_iter0(_iter0)	{
		m_maxiter = 150; m_miniter = 3;
	}
public:
	int GetBundlLineNum(){ return m_vst3d.size(); } //������Ŀ
	void GetSpacePoint3D(SpacePoint3D &spt){ spt.X = m_X; spt.Y = m_Y; spt.Z = m_Z; }
	void GetSpacePoint3D_m0(double &m0) { m0 = m_m0; } //��λȨ����
	void GetSpacePoint3D_Qxx(double *Qxx){ memcpy(Qxx, m_Qxx, 9 * sizeof(double)); }//Э�������
	void GetSpacePoint3D_Dxx(double &Dxx, double &Dyy, double &Dzz){ //����
		Dxx = m_m0 * m_Qxx[0]; Dyy = m_m0 * m_Qxx[4]; Dzz = m_m0 * m_Qxx[8]; }
	void RunintersinglePointForwardCross(){ InterationForwardCross(); }
private:
	void CalRmatrix(Station3D st, double *R);//��ת����
	void InterationForwardCross(); //������ǰ������
	bool DataPrepareWarning(); //׼������
private:
	//������ֵ
	double m_X0; //�������ĳ�ֵ
	double m_Y0;
	double m_Z0;
	double m_iter0; //�������ȣ���XYZһ�£�
	double m_maxiter; //����������
	double m_miniter; //��С��������

	//���ֵ
	double m_X;
	double m_Y;
	double m_Z;

	//��λȨ����, Э�������
	double m_m0;
	double m_Qxx[9];

	vector<Point2D> m_vpt2d;
	vector<Station3D> m_vst3d;
};

#endif