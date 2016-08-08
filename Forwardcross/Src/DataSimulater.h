#ifndef DATASIMULATER_H
#define DATASIMULATER_H
//产生像点模拟数据

#include "DataStruct.h"

double _grnl(double u, double g); //产生一个均值为u,标准差为g的随机数

//正太分布
void _mgrns(double u, double g, /*double *r,*/ int n, double *a);

//均匀分布
void _rnds(double a, double b, double *p, int n);

void UniformDis2D(Point2D pt, double thresold_x, double thresold_y, vector<Point2D> &vpt, int n, RandomSym rsym = Normal_Distribution);

int ModifyPoint2DPair(vector<Point2D> &vpt1, vector<Point2D> &vpt2, int n1, int n2);

void TickPoint2D(vector<Point2D> &vpt, Point2D pt, double threshold);

int SimulatePointPairList(Point2D pl, Point2D pr, vector<Point2D> &pl_list, vector<Point2D> &pr_list, 
	double threshold_px, int threshold_num, RandomSym rsym = Normal_Distribution);

//////////////////////////////////////////////////////////////////////////
//File Operator
void SavePoint2D(const char *file, vector<Point2D> vpt);

void SavePoint2D_I2(const char *file, vector<Point2D> vpt1, vector<Point2D> vpt2);

void SavePoint2D_I1(const char *file, vector<Point2D> vpt, Point2D pt);

void SaveSpacePoint3D(const char *file, vector<SpacePoint3D> vpt3d);

void LoadVirtualPoint(const char *file, vector<WxVirtualPoint> &virtualpoint); //王祥

void LoadVirtulPoint2File(const char *file, vector<SpacePoint3D> &vspt, 
	vector<Point2D> &vptl, vector<Point2D> &vptr);

//////////////////////////////////////////////////////////////////////////

void MulPointSingleForwardCross(vector<Point2D> ptl, vector<Point2D> ptr, 
	Station3D st1, Station3D st2, vector<SpacePoint3D> &vpt3d);

double PointPairLocationSep(Point2D pl, Point2D pr, SpacePoint3D spt, 
	Station3DBLH stl, Station3DBLH str, double threshold_px, int threshold_max, 
	int min_rand_point = 5, double a = WGS84_a, double e2 = WGS84_e2, RandomSym rsym = Normal_Distribution);//一对同名点的SEP

double StatisticPointPairSEP(Point2D pl, Point2D pr, SpacePoint3D spt,
	Station3DBLH stl, Station3DBLH str, double threshold_px, int threshold_max,
	int min_rand_point = 5, int cal_num = 500, double a = WGS84_a, double e2 = WGS84_e2, RandomSym rsym = Normal_Distribution); //基于统计数据计算SEP

//////////////////////////////////////////////////////////////////////////
//统计信息
void StatisticMeanX(vector<SpacePoint3D> spt, double *X);//统计向量的均值

void StatisticVarX(vector<SpacePoint3D> spt, double *var_X); //统计向量的方差

void StatisticMaxX(vector<SpacePoint3D> spt, double *max_X);
void StatisticMinX(vector<SpacePoint3D> spt, double *min_X);

double  CalVectorMax(double *x, int n);
double  CalVectorMin(double *x, int n);

double  DisSpacePoint3D(SpacePoint3D spt1, SpacePoint3D spt2);

double  ApproximateSEP(vector<SpacePoint3D> spt, double sptNum_threshold);//近似计算SEP
double  ApproximateSEP2(vector<SpacePoint3D> spt, SpacePoint3D pt3d, double sptNum_threshold);

double  CalSEP1D(double *x, double *y, double *z, double meanx, double meany, double meanz, int n, double num_threshold);
double* CalSep2D(double *x2D, double *y2D, double *z2D, double *meanx, double *meany, double *meanz, int m, int n, double num_threshold);

double  CalMean1D(double *x, int n);
double* CalMean2D(double *x, int m, int n);
double  CalVar1D(double *x, double mean0, int n);
double* CalVar2D(double *x, double *mean0, int m, int n);
double  DisVector3D(double x1, double y1, double z1, double x2, double y2, double z2);

//////////////////////////////////////////////////////////////////////////
//以下可求解一对同名点范围内的理论最大定位误差

double MaxResidulPointPair(Point2D pl, Point2D pr, Station3D stl, Station3D str, SpacePoint3D spacept3d, 
	double max_threshold_px, double *dX = NULL, FcrossMethod method = ForwardCross_Linear);

Point2D Point2DAddMaxRes(Point2D pt, double threshold);

//////////////////////////////////////////////////////////////////////////
//1. 考虑标称图的姿态测量误差  中误差：20urad 三个测角
//2. 考虑大地坐标误差          中误差 ：B, L 20urad , H 30m
//3. 像点误差                 中误差  ： 60urad(纵横量测值不相关)

double* SymstemResidualForwardCross(Point2D pl, Point2D pr, Station3DBLH st_BLH_L, Station3DBLH st_BLH_R, 
	SpacePoint3D spacept3d, double max_threshold_px, double max_threshold_py, double max_threshold_BL, double max_threshold_H, 
	double max_threshold_phi, double max_threshold_omega, double max_threshold_kappa,/* vector<SpacePoint3D> &vpt3d,*/ int max_point = 5000,
	RandomSym sym =  Normal_Distribution, FcrossMethod method = ForwardCross_Linear, double a = WGS84_a, double e2 = WGS84_e2);

void SimulateStaBLH(Station3DBLH st_BLH, double max_B, double max_L, double max_H,
	double max_p, double max_o, double max_k, int max_num,
	vector<Station3DBLH> &vst_BLH_List, RandomSym rsm = Normal_Distribution);

void StatisticVar2Spt(vector<SpacePoint3D> spt, SpacePoint3D pt3d, double *var_X);



//////////////////////////////////////////////////////////////////////////

void MetoCaloSimulate(vector<Point2D> pl_List, vector<Point2D> pr_List, Station3DBLH st_BLH_L, Station3DBLH st_BLH_R,
	vector<SpacePoint3D> spt_List, double px_g2, double py_g2, double B_g2, double L_g2, double H_g2,
	double p_g2, double o_g2, double k_g2, int max_count = 5000, FILE *fp = NULL,
	RandomSym sym = Normal_Distribution, FcrossMethod method = ForwardCross_Linear, double a = WGS84_a, double e2 = WGS84_e2);

#endif