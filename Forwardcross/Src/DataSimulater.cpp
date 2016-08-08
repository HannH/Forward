#include "DataSimulater.h"

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "SinglePointFcoss.h"
#include "SGCordinate.h"

#pragma warning(disable : 4018)

double _grnl(double u, double g) //产生一个均值为u,方差为g的随机数
{
	double *r, q;
	r = &q;
	srand(GetTickCount()); ////GetTick()（从操作系统到现在的毫秒数）的不同，产生不同的随机数种子值，使得rand函数每次的值都不同
	*r = (double)(rand() % 10);
	
	int i, m;
	double s, w, v, t;
	s = 65536.0; w = 2053.0; v = 13849.0;
	t = 0.0;

	for (i = 1; i <= 12; i++)
	{
		*r = (*r) * w + v; m = (int)(*r / s);
		*r = *r - m * s; t = t + (*r) / s;
	}
	t = u + g * (t - 6.0);

	//Sleep(1);//函数滞留1s
	return t;
}

/************************************************************************/
/* 产生均值为u,方差为g的正太随机分布序列
   u   - 均值
   g   - 方差
   n   - 个数
   a   - 序列
*/
/************************************************************************/
void _mgrns(double u, double g, /*double *r,*/ int n, double *a)
{
	int i = 0;
	double r[12];
	for (i = 0; i < 12; i++)
		r[i] = 0;
	srand(GetTickCount());//GetTick()（从操作系统到现在的毫秒数）的不同，产生不同的随机数种子值，使得rand函数每次的值都不同
	for (i = 0; i < 12; i++)	//注意此处可能要进一步随机
	{
		r[i] = rand() % 10;
	}
	int k, m;
	double s, w, v, t;
	s = 65536.0; w = 2053.0; v = 13849.0;
	for (k = 0; k <= n - 1; k++)
	{
		t = 0.0;
		for (i = 0; i < 12; i++)
		{
			*r = (*r)*w + v;
			m = (int)(*r / s);
			*r = *r - m*s;
			t = t + (*r) / s;
		}
		a[k] = u + g*(t - 6.0);
	}
	Sleep(5);//函数滞留1s
	return ;
}

/************************************************************************/
/* 产生[a,b]之间的均匀分布序列
   a   - 下限
   b   - 上限
   p   - 数组
   n   - 个数*/
/************************************************************************/
void _rnds(double a, double b, double *p, int n)
{
	double *r, q;
	r = &q;
	srand(GetTickCount()); ////GetTick()（从操作系统到现在的毫秒数）的不同，产生不同的随机数种子值，使得rand函数每次的值都不同
	*r = (double)(rand() % 10);
	if (int(*r) % 2 == 1)
		*r += 1; //种子点选奇数

	int i, m;
	double s, u, v;
	s = 65536.0; u = 2053.0; v = 13849.0;
	for (i = 0; i < n; i++)
	{
		*r = u * (*r) + v;
		m = (int)(*r / s);
		*r = *r - m * s;
		p[i] = *r / s;  //此时产生的是[0,1]之间的随机序列

		//进行拉伸
		p[i] = (b - a) * p[i] + a;
	}
	Sleep(5);//函数滞留5ms
	return;
}

/************************************************************************/
/* 在点pt位置附近产生均匀分布的点
   pt            - 原始点位
   threshold_x   - x方向偏离原始点位的距离
   threshold_y   - y方向偏离原始点位的距离
   vpt           - 二维像点数组
   n             - 随机像点数目*/
/************************************************************************/
void UniformDis2D(Point2D pt, double thresold_x, double thresold_y, vector<Point2D> &vpt, int n, RandomSym rsym /*= Normal_Distribution*/)
{
	vpt.clear();

	double a_x, b_x, a_y, b_y;
	a_x = pt.alpha - thresold_x; b_x = pt.alpha + thresold_x;
	a_y = pt.beta - thresold_y ; b_y = pt.beta + thresold_y;

	double *x = new double[n];
	double *y = new double[n];

	switch (rsym){
	case Normal_Distribution:{
								 _mgrns(pt.alpha, 0.25 * thresold_x * thresold_x, n, x); //2dlt <u,  dlt * dlt = 0.25*U*U
								 _mgrns(pt.beta, 0.25 * thresold_y * thresold_y, n, y);
	}
		break;
	case Uniform_Distribution:{
								  _rnds(a_x, b_x, x, n);
								  _rnds(a_y, b_y, y, n);
	}
		break;
	default:{
				_mgrns(pt.alpha, 0.25 * thresold_x * thresold_x, n, x); //2dlt <u,  dlt * dlt = 0.25*U*U
				_mgrns(pt.beta, 0.25 * thresold_y * thresold_y, n, y);
	}
		break;
	}
	
	vpt.resize(n);
	for (int i = 0; i < n; i++)
	{
		vpt[i].imgID = pt.imgID;
		vpt[i].alpha = x[i];
		vpt[i].beta  = y[i];
	}

	delete[] x;
	delete[] y;
}

/************************************************************************/
/* 进行同名点的截断，返回同名点对数
   vpt1 - 左影像序列
   vpt2 - 右影像序列
   n1   - 左影像点数
   n2   - 右影像像点数*/
/************************************************************************/
int ModifyPoint2DPair(vector<Point2D> &vpt1, vector<Point2D> &vpt2, int n1, int n2)
{
	int n = min(n1, n2);
	vector<Point2D> pt;
	pt.resize(n);

	if (n1 < n2)
	{//将n2 截断
		for (int i = 0; i < n; i++)
		{
			pt[i].imgID = vpt2[i].imgID;
			pt[i].alpha = vpt2[i].alpha;
			pt[i].beta = vpt2[i].beta;
		}
		vpt2.clear();
		vpt2.resize(n);
		vpt2.assign(pt.begin(), pt.end());
		//pt.clear();
	}
	else{//将n1截断
		for (int i = 0; i < n; i++)
		{
			pt[i].imgID = vpt1[i].imgID;
			pt[i].alpha = vpt1[i].alpha;
			pt[i].beta = vpt1[i].beta;
		}

		vpt1.clear();
		vpt1.resize(n);
		vpt1.assign(pt.begin(), pt.end());
		//pt.clear();
	}

	//pt.clear();
	return n;
}

/************************************************************************/
/* 剔除中心像点误差限之外的点  
  vpt       -  点序列         ||dc|| <= threshold
  pt        -  中心像点
  threshold -  阈值（弧度）  tan(a)*tan(a) + tan(b)*tan(b) = tan(c)* tan(c)_*/
/************************************************************************/
void TickPoint2D(vector<Point2D> &vpt, Point2D pt, double threshold)
{
	vector<Point2D> pt_new;
	int num = 0;

	//double a, b, a0, b0;
	//double tc2, tc2_new; //tan(c) * tan(c)
	//double cc2, cc2_new;
	//double dlt_c;
	//double l1, l2, l3; //l1 * l1 + l2 * l2 - 2 * l1 * l2 * cos(dlt_c)= l3 * l3

	//a0 = pt.alpha; b0 = pt.beta;
	//int num = 0;
	//for (int i = 0; i < vpt.size(); i++)
	//{
	//	tc2 = tan(a0) * tan(a0) + tan(b0) * tan(b0);
	//	cc2 = 1 / (1 + tc2);

	//	a = vpt[i].alpha; b = vpt[i].beta;
	//	tc2_new = tan(a) * tan(a) + tan(b) * tan(b);
	//	cc2_new = 1 / (1 + tc2_new);
	//	
	//	l3 =  (tan(a) - tan(a0)) * (tan(a) - tan(a0)) + (tan(b) - tan(b0)) * (tan(b) - tan(b0));
	//	l1 = 1 / cc2;
	//	l2 = 1 / cc2_new;
	//	dlt_c = 0.5 * (l1 + l2 - l3) / sqrt(l1 * l2);
	//	dlt_c = acos(dlt_c);

	//	if (fabs(dlt_c) <= threshold)
	//	{
	//		pt_new.push_back(vpt[i]); num++;
	//	}
	//}

	for (int i = 0; i < vpt.size(); i++)
	{
		if (fabs(vpt[i].alpha - pt.alpha) <= threshold && fabs(vpt[i].beta - pt.beta) <= threshold) {
			pt_new.push_back(vpt[i]); num++;
		}	
	}

	vpt.clear();
	vpt.resize(num);
	vpt.assign(pt_new.begin(), pt_new.end());

	pt_new.clear();
}

/************************************************************************/
/* 同名点附近产生随机点
   pl, pr             - 同名点对
   pl_list, pr_list   - 产生的序列点对
   threshold_px       - 像素阈值
   threshold_num      - 点数阈值*/
/************************************************************************/
int SimulatePointPairList(Point2D pl, Point2D pr, vector<Point2D> &pl_list, vector<Point2D> &pr_list,
	double threshold_px, int threshold_num, RandomSym rsym /*= Normal_Distribution*/)
{
	UniformDis2D(pl, threshold_px, threshold_px, pl_list, threshold_num, rsym); 
	//TickPoint2D(pl_list, pl, threshold_px);
	UniformDis2D(pr, threshold_px, threshold_px, pr_list, threshold_num, rsym); 
	//TickPoint2D(pr_list, pr, threshold_px);

	int n1 = pl_list.size();
	int n2 = pr_list.size();

	return threshold_num;
	//return ModifyPoint2DPair(pl_list, pr_list, n1, n2);
}

//保存单点文件
void SavePoint2D(const char *file, vector<Point2D> vpt)
{
	FILE *fp = fopen(file, "w");
	if (!fp)
	{
		printf("保存像点文件有误!\n");
		exit(0);
	}
	fprintf(fp, "$$  Point Number\n");
	fprintf(fp, "$$  lx   ly\n");
	fprintf(fp, "%d\n", vpt.size());
	for (int i = 0; i < vpt.size(); i++)
	{
		fprintf(fp, "%.10e   %.10e\n", vpt[i].alpha, vpt[i].beta);
	}
	fclose(fp);
}

//保存同名点文件
void SavePoint2D_I2(const char *file, vector<Point2D> vpt1, vector<Point2D> vpt2)
{
	int num = min(vpt1.size(), vpt2.size());

	FILE *fp = fopen(file, "w");
	if (!fp)
	{
		printf("保存像点文件有误!\n");
		exit(0);
	}
	fprintf(fp, "$$  Point Pair Number\n");
	fprintf(fp, "$$  lx   ly      rx   ry\n");
	fprintf(fp, "%d\n", num);
	for (int i = 0; i < num; i++)
	{
		fprintf(fp, "%.10e   %.10e      %.10e   %.10e\n", vpt1[i].alpha, vpt1[i].beta, vpt2[i].alpha, vpt2[i].beta);
	}
	fclose(fp);
}

//保存中心点附近的像点文件，以及残差
void SavePoint2D_I1(const char *file, vector<Point2D> vpt, Point2D pt)
{
	int num = vpt.size();

	FILE *fp = fopen(file, "w");
	if (!fp)
	{
		printf("保存像点文件有误!\n");
		exit(0);
	}
	fprintf(fp, "$$  Point Number\n");
	fprintf(fp, "$$  Orignal Point\n");
	fprintf(fp, "$$  lx   ly      rx   ry\n");
	fprintf(fp, "%d\n", num);
	fprintf(fp, "%.10e  %.10e\n", pt.alpha, pt.beta);
	for (int i = 0; i < num; i++)
	{
		fprintf(fp, "%.10e   %.10e      %.10e   %.10e\n", 
			vpt[i].alpha, vpt[i].beta, 
			vpt[i].alpha - pt.alpha, vpt[i].beta - pt.beta);
	}
	fclose(fp);
}

//保存物方点文件
void SaveSpacePoint3D(const char *file, vector<SpacePoint3D> vpt3d)
{
	FILE *fp = fopen(file, "w");
	if (!fp)
	{
		printf("保存像点文件有误!\n");
		exit(0);
	}
	fprintf(fp, "$$  Space Point Number\n");
	fprintf(fp, "$$  ID   X    Y    Z\n");
	fprintf(fp, "%d\n", vpt3d.size());
	for (int i = 0; i < vpt3d.size(); i++)
	{
		vpt3d[i].spaceID = i;
		fprintf(fp, "%d   %.6lf   %.6lf    %.6lf\n", vpt3d[i].spaceID, vpt3d[i].X, vpt3d[i].Y, vpt3d[i].Z);
	}
	fclose(fp);
}

//加载模拟轨道点
void LoadVirtualPoint(const char *file, vector<WxVirtualPoint> &virtualpoint)
{
	FILE *fp = fopen(file, "r");
	if (!fp)
	{
		printf("读取虚拟文件有误!\n"); exit(0);
	}
	char buf[512];
	fgets(buf, 512, fp);

	int num = 0;
	WxVirtualPoint vp;
	while (1) {
		if (feof(fp)) {
			break;
		}
		fgets(buf, 512, fp);
		int ID;
		double time, X, Y, Z, lx, ly, rx, ry;
		num = sscanf(buf, "%d%lf%lf%lf%lf%lf%lf%lf%lf", &ID, &time, &X, &Y, &Z, &lx, &ly, &rx, &ry);
		if (num < 9) continue;
		vp.ID = ID; vp.time = time; 
		vp.pt3d.X = X; vp.pt3d.Y = Y; vp.pt3d.Z = Z;
		vp.pl.alpha = lx; vp.pl.beta = ly;
		vp.pr.alpha = rx; vp.pr.beta = ry;
		virtualpoint.push_back(vp);
	}
	fclose(fp);
}

void LoadVirtulPoint2File(const char *file, vector<SpacePoint3D> &vspt,
	vector<Point2D> &vptl, vector<Point2D> &vptr)
{
	FILE *fp = fopen(file, "r");
	if (!fp)
	{
		printf("读取虚拟文件有误!\n"); exit(0);
	}
	char buf[512];
	fgets(buf, 512, fp);

	int num = 0;
	SpacePoint3D spt;
	Point2D pl, pr;
	while (1) {
		if (feof(fp)) {
			break;
		}
		fgets(buf, 512, fp);
		int ID;
		double time, X, Y, Z, lx, ly, rx, ry;
		num = sscanf(buf, "%d%lf%lf%lf%lf%lf%lf%lf%lf", &ID, &time, &X, &Y, &Z, &lx, &ly, &rx, &ry);
		if (num < 9) continue;
		//vp.ID = ID; vp.time = time;
		spt.X = X; spt.Y = Y; spt.Z = Z;
		pl.alpha = lx; pl.beta = ly; pl.imgID = 0;
		pr.alpha = rx; pr.beta = ry; pr.imgID = 1;
		vspt.push_back(spt);
		vptl.push_back(pl); vptr.push_back(pr);
	}
	fclose(fp);

}

//////////////////////////////////////////////////////////////////////////

//通过多组数据进行前方交会，结果保存于vpt3d中
void MulPointSingleForwardCross(vector<Point2D> ptl, vector<Point2D> ptr,
	Station3D st1, Station3D st2, vector<SpacePoint3D> &vpt3d)
{
	vpt3d.clear();

	if (ptl.size() != ptr.size())
	{
		printf("同名点数据不对应! n1 = %d, n2 = %d\n", ptl.size(), ptr.size());
		system("pause");
		exit(0);
	}

	for (int i = 0; i < ptl.size(); i++)
	{
		vector<Point2D> pointpair;
		vector<Station3D> stationpair;
		
		pointpair.push_back(ptl[i]); pointpair.push_back(ptr[i]);
		stationpair.push_back(st1) ; stationpair.push_back(st2);

		CsinglePointForwardCross sf(pointpair, stationpair);
		sf.RunSinglePointForwardCross();
		SpacePoint3D pt3d;
		sf.GetSpacePoint3D(pt3d);
		vpt3d.push_back(pt3d);
	}
}


//////////////////////////////////////////////////////////////////////////
//统计信息
//统计向量的均值
void StatisticMeanX(vector<SpacePoint3D> spt, double *X)
{
	X[0] = X[1] = X[2] = 0.0;
	int n = spt.size();
	for (int i = 0; i < n; i++)
	{
		X[0] += spt[i].X; X[1] += spt[i].Y; X[2] = spt[i].Z;
	}
	X[0] /= n; X[1] /= n; X[2] /= n;
}

//统计向量的方差
void StatisticVarX(vector<SpacePoint3D> spt, double *var_X)
{
	if (spt.size() <= 0) {
		var_X[0] = var_X[1] = var_X[2] = 0.0;
		return;
	}

	double meax_X[3];
	StatisticMeanX(spt, meax_X);

	int n = spt.size();
	var_X[0] = var_X[1] = var_X[2] = 0.0;
	for (int i = 0; i < n; i++)
	{
		var_X[0] += (spt[i].X - meax_X[0]) * (spt[i].X - meax_X[0]);
		var_X[1] += (spt[i].Y - meax_X[1]) * (spt[i].Y - meax_X[1]);
		var_X[2] += (spt[i].Z - meax_X[2]) * (spt[i].X - meax_X[2]);
	}

	var_X[0] /= (n - 1);  //为什么（n-1）
	var_X[1] /= (n - 1);
	var_X[2] /= (n - 1);
}

void StatisticMaxX(vector<SpacePoint3D> spt, double *max_X)
{
	max_X[0] = max_X[1] = max_X[2] = DOUBLE_MAX;
	for (int i = 0; i < spt.size(); i++)
	{
		if (spt[i].X > max_X[0]) max_X[0] = spt[i].X;
		if (spt[i].Y > max_X[1]) max_X[1] = spt[i].Y;
		if (spt[i].Z > max_X[2]) max_X[2] = spt[i].Z;
	}
}

void StatisticMinX(vector<SpacePoint3D> spt, double *min_X)
{
	min_X[0] = min_X[1] = min_X[2] = DOUBLE_MIN;
	for (int i = 0; i < spt.size(); i++)
	{
		if (spt[i].X < min_X[0]) min_X[0] = spt[i].X;
		if (spt[i].Y < min_X[1]) min_X[1] = spt[i].Y;
		if (spt[i].Z < min_X[2]) min_X[2] = spt[i].Z;
	}
}

double DisSpacePoint3D(SpacePoint3D spt1, SpacePoint3D spt2)
{
	double r;
	r = (spt1.X - spt2.X) * (spt1.X - spt2.X) + (spt1.Y - spt2.Y) * (spt1.Y - spt2.Y) + (spt1.Z - spt2.Z) * (spt1.Z - spt2.Z);
	return sqrt(r);
}

double  CalVectorMax(double *x, int n)
{
	double amax = DOUBLE_MAX;
	for (int i = 0; i < n; i++) {
		if (x[i] > amax) {
			amax = x[i];
		}
	}
	return amax;
}
double  CalVectorMin(double *x, int n)
{
	double amin = DOUBLE_MIN;
	for (int i = 0; i < n; i++) {
		if (x[i] < amin) {
			amin = x[i];
		}
	}
	return amin;
}
double  CalMean1D(double *x, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i];
	return (sum / n);
}
double* CalMean2D(double *x, int m, int n)//行向量的均值, m行，n列
{
	double *sum = new double[m];
	memset(sum, 0, m * sizeof(double));

	double *x_l;
	x_l = x;
	for (int i = 0; i < m; i ++)	{
		for (int j = 0; j < n; j++) {
			sum[i] += x_l[j];
		}
		sum[i] /= n;
     	x_l += n;
	}
	return sum;
}
double  CalVar1D(double *x, double mean0, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)	{
		sum += (x[i] - mean0) * (x[i] - mean0);
	}
	sum = sqrt(sum / n);
	return sum;
}
double* CalVar2D(double *x, double *mean0, int m, int n)//行向量的方差, m行，n列
{
	double *Var = new double[m];
	memset(Var, 0, m * sizeof(double));

	double *x_l;
	x_l = x;
	for (int i = 0; i < m; i++)	{
		for (int j = 0; j < n; j++)	{
			Var[i] += (x_l[j] - mean0[i]) * (x_l[j] - mean0[i]);
		}
		Var[i] = sqrt(Var[i] / n);
		x_l += n;
	}
	return Var;
}

double DisVector3D(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double r = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
	return (sqrt(r));
}


//////////////////////////////////////////////////////////////////////////
//SEP
//近似计算SEP
double ApproximateSEP(vector<SpacePoint3D> spt, double sptNum_threshold)
{
	if (spt.size() < sptNum_threshold) {
		printf("数据量太少，不利于统计分析！ %d\n", sptNum_threshold);
		return 0;
	}
	double mean_X[3];
	StatisticMeanX(spt, mean_X); //重心坐标

	SpacePoint3D spt0(mean_X[0], mean_X[1], mean_X[2]);

	double max_X[3], min_X[3];
	StatisticMaxX(spt, max_X);
	StatisticMinX(spt, min_X);

	double max_R = max((max_X[0] - min_X[0]), (max_X[1] - min_X[1]));
	max_R = 0.5 * max((max_X[2] - min_X[2]), max_R);

	//printf("maxR : %lf\n", max_R);

	double sep_p = 1.0; //SEP的概率
	double sep_threshold = 0.5; //SEP的阈值
	int sep_num = 0;

	int n = spt.size();
	double *r = new double[n];

	//printf("到重心的距离: %lf %lf %lf\n", mean_X[0], mean_X[1], mean_X[2]);
	for (int i = 0; i < n; i++) {
		r[i] = DisSpacePoint3D(spt0, spt[i]);
	}

	double min_R = 0;
	double R = max_R; //搜索半径
	double R_threshold = 0.005;//5米的阈值
	do 
	{
		for (int i = 0; i < n; i++)	{
			if (r[i] < R) sep_num++;
		}

		sep_p = (double)sep_num / n; //计算概率
		if (sep_p >= sep_threshold) {//半径可进一步缩小
			max_R = R; //下次迭代的上限
			R = 0.5 * (R + min_R); //搜索半径减小一半
		}
		if (sep_p < sep_threshold) {//半径可扩大
			min_R = R; //下次迭代的下限
			R = 0.5 * (max_R + R);
		}
	} while (fabs(R - min_R) <= R_threshold && fabs(R - max_R) <= R_threshold);

	//printf("SEP: %lf\n", R);
	delete[] r;
	return R;
}

double ApproximateSEP2(vector<SpacePoint3D> spt, SpacePoint3D pt3d, double sptNum_threshold)
{
	if (spt.size() < sptNum_threshold) {
		printf("数据量太少，不利于统计分析！ %d\n", sptNum_threshold);
		return 0;
	}
	double mean_X[3];
	StatisticMeanX(spt, mean_X); //重心坐标

	printf("重心: %lf %lf %lf\n", mean_X[0], mean_X[1], mean_X[2]);
	SpacePoint3D spt0(mean_X[0], mean_X[1], mean_X[2]);

	double max_X[3], min_X[3];
	StatisticMaxX(spt, max_X);
	StatisticMinX(spt, min_X);

	double max_R = max((max_X[0] - min_X[0]), (max_X[1] - min_X[1]));
	max_R = 0.5 * max((max_X[2] - min_X[2]), max_R);

	printf("maxR : %lf\n", max_R);
	double sep_p = 1.0; //SEP的概率
	double sep_threshold = 0.5; //SEP的阈值
	int sep_num = 0;

	int n = spt.size();
	double *r = new double[n];

	for (int i = 0; i < n; i++) {
		r[i] = DisSpacePoint3D(pt3d, spt[i]);
	}

	double min_R = 0;
	double R = max_R; //搜索半径
	double R_threshold = 1e-6;//3米的阈值
	do
	{
		for (int i = 0; i < n; i++)	{
			if (r[i] <= R) sep_num++;
		}

		sep_p = (double)sep_num / n; //计算概率
		if (sep_p >= sep_threshold) {//半径可进一步缩小
			max_R = R; //下次迭代的上限
			R = 0.5 * (R + min_R); //搜索半径减小一半
		}
		if (sep_p < sep_threshold) {//半径可扩大
			min_R = R; //下次迭代的下限
			R = 0.5 * (max_R + R);
		}
	} while (fabs(R - min_R) <= R_threshold && fabs(R - max_R) <= R_threshold);

	//printf("SEP: %lf\n", R);
	delete[] r;
	return R;
}

//一个点位的SEP
double CalSEP1D(double *x, double *y, double *z, double meanx, double meany, double meanz, int n, double num_threshold)
{
	if (n < num_threshold) {
		printf("数据量太少，不利于统计分析！ %d\n", n);
		return 0;
	}

	double min_x = CalVectorMin(x, n);
	double min_y = CalVectorMin(y, n);
	double min_z = CalVectorMin(z, n);
	double max_x = CalVectorMax(x, n);
	double max_y = CalVectorMax(y, n);
	double max_z = CalVectorMax(z, n);

	double max_R = max((max_x - min_x), (max_y - min_y));
	max_R = max((max_z - min_z), max_R); //选择直径作为迭代阈值

	double *r = new double[n];
	for (int i = 0; i < n; i++) {
		r[i] = DisVector3D(x[i], y[i], z[i], meanx, meany, meanz);
	}

	double sep_p = 1.0; //SEP的概率
	double sep_threshold = 0.5; //SEP的阈值
	int sep_num = 0;
	double min_R = 0;
	double R = max_R; //搜索半径
	double R_threshold = 1e-6;//1米的阈值
	do
	{
		for (int i = 0; i < n; i++)	{
			if (r[i] <= R) sep_num++;
		}

		sep_p = (double)sep_num / n; //计算概率
		if (sep_p >= sep_threshold) {//半径可进一步缩小
			max_R = R; //下次迭代的上限
			R = 0.5 * (R + min_R); //搜索半径减小一半
		}
		if (sep_p < sep_threshold) {//半径可扩大
			min_R = R; //下次迭代的下限
			R = 0.5 * (max_R + R);
		}
	} while (fabs(R - min_R) <= R_threshold && fabs(R - max_R) <= R_threshold);

	delete[] r;
	return R;
}

//多个点位的SEP
double* CalSep2D(double *x2D, double *y2D, double *z2D, double *meanx, double *meany, double *meanz, int m, int n, double num_threshold)
{
	if (n < num_threshold) {
		printf("数据量太少，不利于统计分析！ %d\n", n);
		return 0;
	}

	double *SEP = new double[m];
	double *x, *y, *z;
	x = x2D; y = y2D; z = z2D; 
	for (int i = 0; i < m; i++)
	{
		SEP[i] = CalSEP1D(x, y, z, meanx[i], meany[i], meanz[i], n, num_threshold);
		x += n; y += n; z += n; 
	}
	return SEP;
}


/************************************************************************/
/* pl， pr        -    真实像点
   stl, str       -    左右摄站真实的位置
   spt            -    真实的空间坐标
   threshold_px   -    像方的误差阈值
   threshold_max  -    随机点的最大值
   min_rand_point -    随机点的最小值
   a, e2          -    椭球参数*/
/************************************************************************/
double PointPairLocationSep(Point2D pl, Point2D pr, SpacePoint3D spt,
	Station3DBLH stl, Station3DBLH str, double threshold_px, int threshold_max, 
	int min_rand_point/* = 5*/, double a/* = WGS84_a*/, double e2 /*= WGS84_e2*/, 
	RandomSym rsym/* = Normal_Distribution*/)
{
	vector<Point2D> vpl_List, vpr_List;
	double percent = 0.0;

	int ncount = 0; //循环次数
	do 
	{//最多循环10次
		int sim_num = SimulatePointPairList(pl, pr, vpl_List, vpr_List, threshold_px, threshold_max, rsym);
		percent = (double)sim_num / threshold_max;
		//printf("个数 = %d, 比例 = %lf\n", sim_num, percent);
		if (percent >= 0.5) break;
		ncount++;
	} while (ncount < 10);
	
	Station3D Sta_L, Sta_R;
	StaBLH2XYZ(stl, Sta_L, a, e2);
	StaBLH2XYZ(str, Sta_R, a, e2);

	vector<SpacePoint3D> vpt3d;	
	MulPointSingleForwardCross(vpl_List, vpr_List, Sta_L, Sta_R, vpt3d);

	SaveSpacePoint3D("space3d.txt", vpt3d);

	double sep = ApproximateSEP2(vpt3d, spt, min_rand_point);

	printf("sep = %lf\n", sep);

	vpl_List.clear();
	vpr_List.clear();
	vpt3d.clear();

	return sep;
}

/************************************************************************/
/*    基于统计数据计算SEP 
pl， pr        -    真实像点
stl, str       -    左右摄站真实的位置
spt            -    真实的空间坐标
threshold_px   -    像方的误差阈值
threshold_max  -    随机点的最大值
min_rand_point -    随机点的最小值
a, e2          -    椭球参数* /
/************************************************************************/
double StatisticPointPairSEP(Point2D pl, Point2D pr, SpacePoint3D spt,
	Station3DBLH stl, Station3DBLH str, double threshold_px, int threshold_max,
	int min_rand_point/* = 5*/, int cal_num /*= 500*/, double a /*= WGS84_a*/, double e2 /*= WGS84_e2*/,
	RandomSym rsym/* = Normal_Distribution*/)

{
	double sep[500] = { 0.0 };
	double sum = 0.0;
	for (int i = 0; i < cal_num; i++) {
		sep[i] = PointPairLocationSep(pl, pr, spt, stl, str, threshold_px, threshold_max, min_rand_point, a, e2, rsym);
		sum += sep[i];
	}
	
	sum /= cal_num;
	return sum;
}


//////////////////////////////////////////////////////////////////////////
/************************************************************************/
/* 计算单点理论定位最大误差
   pl, pr                        -- 同名点
   stl, str,                     -- 摄站点
   pt3d                          -- 真实空间点
   double                        -- 检查最大误差
   dX                            -- 三个方向误差*/
/************************************************************************/
double MaxResidulPointPair(Point2D pl, Point2D pr, Station3D stl, Station3D str, SpacePoint3D spacept3d,
	double max_threshold_px, double *dX/* = NULL*/, FcrossMethod method /*= ForwardCross_Linear*/)
{
	max_threshold_px = fabs(max_threshold_px);
	Point2D pl_new, pr_new;
	pl_new = Point2DAddMaxRes(pl, max_threshold_px);
	pr_new = Point2DAddMaxRes(pr, max_threshold_px);

	vector<Point2D> pointpair;
	vector<Station3D> stationpair;

	pointpair.push_back(pl_new); pointpair.push_back(pr_new);
	stationpair.push_back(stl); stationpair.push_back(str);
	SpacePoint3D pt3d;
	switch (method){
	case ForwardCross_Linear: {
								 CsinglePointForwardCross sf(pointpair, stationpair);
								 sf.RunSinglePointForwardCross();
								 sf.GetSpacePoint3D(pt3d);

	}
		break;
	case ForwardCross_Inter: {
								 CintersinglePointForwardCross sisf(pointpair, stationpair,
									 spacept3d.X+10, spacept3d.Y+10, spacept3d.Z+10, 0.001);//初值给10km
								 sisf.RunintersinglePointForwardCross();
								 sisf.GetSpacePoint3D(pt3d);
	}
		break;
	default: break;
	}
	
	if (dX != NULL) {
		dX[0] = spacept3d.X - pt3d.X;
		dX[1] = spacept3d.Y - pt3d.Y;
		dX[2] = spacept3d.Z - pt3d.Z;
	}

	double r = (spacept3d.X - pt3d.X) * (spacept3d.X - pt3d.X) 
		+ (spacept3d.Y - pt3d.Y) * (spacept3d.Y - pt3d.Y) 
		+ (spacept3d.Z - pt3d.Z) * (spacept3d.Z - pt3d.Z);

	return sqrt(r);
}

Point2D Point2DAddMaxRes(Point2D pt, double threshold)
{
	threshold = fabs(threshold);
	Point2D pt_new(pt.imgID, 0, 0);
	pt_new.imgID = pt.imgID; 

	double ta2 = tan(pt.alpha) * tan(pt.alpha);
	double tb2 = tan(pt.beta)  * tan(pt.beta);
	double /*c, cnew,*/ a_new, b_new;

	/*c = atan(sqrt(ta2 + tb2));
	cnew = c + threshold;

	a_new = tan(cnew) * tan(cnew) * ta2 / (ta2 + tb2);
	b_new = tan(cnew) * tan(cnew) * tb2 / (ta2 + ta2);*/

	if (pt.alpha >= 0.0 && pt.beta > 0.0) { //第一象限	
		/*a_new = atan(sqrt(a_new));
		b_new = atan(sqrt(b_new));*/
		a_new = pt.alpha + threshold; b_new = pt.beta + threshold;
	}
	if (pt.alpha < 0.0  && pt.beta >= 0.0) { //第二象限
		/*a_new = -atan(sqrt(a_new));
		b_new = atan(sqrt(b_new));*/
		a_new = pt.alpha - threshold; b_new = pt.beta + threshold;
	}
	if (pt.alpha <= 0.0 && pt.beta <= 0.0) { //第三象限
		/*a_new = -atan(sqrt(a_new));
		b_new = -atan(sqrt(b_new));*/
		a_new = pt.alpha - threshold; b_new = pt.beta - threshold;
	}
	if (pt.alpha > 0.0  && pt.beta < 0.0) { //第四象限
		/*a_new = atan(sqrt(a_new));
		b_new = -atan(sqrt(b_new));*/
		a_new = pt.alpha + threshold; b_new = pt.beta - threshold;
	}

	pt_new.alpha = a_new; pt_new.beta = b_new;

	return pt_new;
}


//////////////////////////////////////////////////////////////////////////

void SimulateStaBLH(Station3DBLH st_BLH, double max_B, double max_L, double max_H,
	double max_p, double max_o, double max_k, int max_num, 
	vector<Station3DBLH> &vst_BLH_List, RandomSym rsm /*= Normal_Distribution*/)
{
	double B = st_BLH.B, L = st_BLH.L, H = st_BLH.H;
	double p = st_BLH.phi, o = st_BLH.omega, k = st_BLH.kappa;

	double *B_new = new double[max_num];
	double *L_new = new double[max_num];
	double *H_new = new double[max_num];
	double *p_new = new double[max_num];
	double *o_new = new double[max_num];
	double *k_new = new double[max_num];

	vst_BLH_List.resize(max_num);

	switch (rsm) {
	case Normal_Distribution:{
								 _mgrns(B, 0.25 * max_B * max_B, max_num, B_new); 
								 _mgrns(L, 0.25 * max_L * max_L, max_num, L_new);
								 _mgrns(H, 0.25 * max_H * max_H, max_num, H_new);
								 _mgrns(p, 0.25 * max_p * max_p, max_num, p_new);
								 _mgrns(o, 0.25 * max_o * max_o, max_num, o_new);
								 _mgrns(k, 0.25 * max_k * max_o, max_num, k_new);
								 for (int i = 0; i < max_num; i++)
								 {
									 vst_BLH_List[i].imgID = st_BLH.imgID; vst_BLH_List[i].rotate = st_BLH.rotate;
									 vst_BLH_List[i].B = B_new[i]; vst_BLH_List[i].L = L_new[i]; vst_BLH_List[i].H = H_new[i];
									 vst_BLH_List[i].phi = p_new[i]; vst_BLH_List[i].omega = o_new[i]; vst_BLH_List[i].kappa = k_new[i];
								 }
	}
		break;
	case Uniform_Distribution: {
								   _rnds(B - max_B, B + max_B, B_new, max_num);
								   _rnds(L - max_L, L + max_L, L_new, max_num);
								   _rnds(H - max_H, H + max_H, H_new, max_num);
								   _rnds(p - max_p, p + max_p, p_new, max_num);
								   _rnds(o - max_o, o + max_o, o_new, max_num);
								   _rnds(k - max_k, k + max_k, k_new, max_num);
								   for (int i = 0; i < max_num; i++)
								   {
									   vst_BLH_List[i].imgID = st_BLH.imgID; vst_BLH_List[i].rotate = st_BLH.rotate;
									   vst_BLH_List[i].B = B_new[i]; vst_BLH_List[i].L = L_new[i]; vst_BLH_List[i].H = H_new[i];
									   vst_BLH_List[i].phi = p_new[i]; vst_BLH_List[i].omega = o_new[i]; vst_BLH_List[i].kappa = k_new[i];
								   }
	}
		break;
	default: break;
	}

	delete[] B_new;
	delete[] L_new;
	delete[] H_new;
	delete[] p_new;
	delete[] o_new;
	delete[] k_new;
}

//蒙特卡洛模拟
double* SymstemResidualForwardCross(Point2D pl, Point2D pr, Station3DBLH st_BLH_L, Station3DBLH st_BLH_R,
	SpacePoint3D spacept3d, double max_threshold_px, double max_threshold_py, double max_threshold_BL, double max_threshold_H,
	double max_threshold_phi, double max_threshold_omega, double max_threshold_kappa, /*vector<SpacePoint3D> &vpt3d,*/ int max_point/* = 5000*/,
	RandomSym sym /*= Normal_Distribution*/, FcrossMethod method /*= ForwardCross_Linear*/, double a/* = WGS84_a*/, double e2/* = WGS84_e2*/)
{
	vector<SpacePoint3D> vpt3d;
	//vpt3d.clear();

	vector<Point2D> pl_List, pr_List;
	vector<Station3DBLH> st_BLH_L_List, st_BLH_R_List;
	UniformDis2D(pl, max_threshold_px, max_threshold_py, pl_List, max_point);
	UniformDis2D(pr, max_threshold_px, max_threshold_py, pr_List, max_point);

	SimulateStaBLH(st_BLH_L, max_threshold_BL, max_threshold_BL, max_threshold_H,
		max_threshold_phi, max_threshold_omega, max_threshold_kappa, max_point, st_BLH_L_List, Normal_Distribution);
	SimulateStaBLH(st_BLH_R, max_threshold_BL, max_threshold_BL, max_threshold_H,
		max_threshold_phi, max_threshold_omega, max_threshold_kappa, max_point, st_BLH_R_List, Normal_Distribution);

	//开始前方交会
	SpacePoint3D spacept3d_new;
	Station3D st_L, st_R;
	for (int i = 0; i < max_point; i++)
	{
		vector<Point2D> pointpair;
		vector<Station3D> stationpair;
		StaBLH2XYZ(st_BLH_L_List[i], st_L, a, e2);
		StaBLH2XYZ(st_BLH_R_List[i], st_R, a, e2);

		pointpair.push_back(pl_List[i]); pointpair.push_back(pr_List[i]);
		stationpair.push_back(st_L); stationpair.push_back(st_R);
		switch (method) {
		case ForwardCross_Linear: {
									  CsinglePointForwardCross sf(pointpair, stationpair);
									  sf.RunSinglePointForwardCross();
									  sf.GetSpacePoint3D(spacept3d_new);
									  vpt3d.push_back(spacept3d_new);
		}
			break;
		case ForwardCross_Inter: {
									 CintersinglePointForwardCross sisf(pointpair, stationpair,
										 spacept3d.X + 10, spacept3d.Y + 10, spacept3d.Z + 10, 0.001);//初值给10km
									 sisf.RunintersinglePointForwardCross();
									 sisf.GetSpacePoint3D(spacept3d_new);
									 vpt3d.push_back(spacept3d_new);
		}
			break;
		default: break;
		}
	}

	//计算sep
	double *Var = new double[3];
	memset(Var, 0, 3 * sizeof(double));
	StatisticVar2Spt(vpt3d, spacept3d, Var);
	/*for (int i = 0; i < 3; i++)
	{
	printf("Var[%d] = %lf\n", i, Var[i]);
	}*/
	
	double sep = ApproximateSEP2(vpt3d, spacept3d, 100);
	return Var;
}

void StatisticVar2Spt(vector<SpacePoint3D> spt, SpacePoint3D pt3d, double *var_X)
{
	if (spt.size() <= 0) {
		var_X[0] = var_X[1] = var_X[2] = 0.0;
		return;
	}

	int n = spt.size();
	var_X[0] = var_X[1] = var_X[2] = 0.0;
	for (int i = 0; i < n; i++)
	{
		var_X[0] += (spt[i].X - pt3d.X) * (spt[i].X - pt3d.X);
		var_X[1] += (spt[i].Y - pt3d.Y) * (spt[i].Y - pt3d.Y);
		var_X[2] += (spt[i].Z - pt3d.Z) * (spt[i].Z - pt3d.Z);
	}

	var_X[0] = sqrt(var_X[0] / (n ));  //为什么（n-1）
	var_X[1] = sqrt(var_X[1] / (n));
	var_X[2] = sqrt(var_X[2] / (n));	
}


//////////////////////////////////////////////////////////////////////////
//第二种方式模拟

void MetoCaloSimulate(vector<Point2D> pl_List, vector<Point2D> pr_List, Station3DBLH st_BLH_L, Station3DBLH st_BLH_R, 
	vector<SpacePoint3D> spt_List, double px_g2, double py_g2, double B_g2, double L_g2, double H_g2, 
	double p_g2, double o_g2, double k_g2, int max_count/* = 5000*/, FILE *fp /*= NULL*/,
	RandomSym sym /*= Normal_Distribution*/, FcrossMethod method /*= ForwardCross_Linear*/, double a/* = WGS84_a*/, double e2/* = WGS84_e2*/)
{
	int point_num = pl_List.size();
	//point_num = 1;

	Station3D stl_L, str_R;
	SpacePoint3D pt_new;

	double tmp_a, tmp_b;
	double tmp_B, tmp_L, tmp_H;
	double tmp_p, tmp_o, tmp_k;

	double *X_List = new double[point_num * max_count]; memset(X_List, 0, point_num * max_count * sizeof(double));
	double *Y_List = new double[point_num * max_count]; memset(Y_List, 0, point_num * max_count * sizeof(double));
	double *Z_List = new double[point_num * max_count]; memset(Z_List, 0, point_num * max_count * sizeof(double));	

	double *MeanX = new double[point_num]; memset(MeanX, 0, point_num * sizeof(double));
	double *MeanY = new double[point_num]; memset(MeanY, 0, point_num * sizeof(double));
	double *MeanZ = new double[point_num]; memset(MeanZ, 0, point_num * sizeof(double));

	for (int i = 0; i < point_num; i++) {
		MeanX[i] = spt_List[i].X; MeanY[i] = spt_List[i].Y; MeanZ[i] = spt_List[i].Z;
	}

	for (int k = 0; k < max_count; k++)
	{
		/*	do { tmp_B = _grnl(st_BLH_L.B, B_g2); } while (tmp_B > 4 * B_g2);
			do { tmp_L = _grnl(st_BLH_L.L, L_g2); } while (tmp_L > 4 * L_g2);
			do { tmp_H = _grnl(st_BLH_L.H, H_g2); } while (tmp_H > 4 * H_g2);
			do { tmp_p = _grnl(st_BLH_L.phi, p_g2); } while (tmp_p > 4 * p_g2);
			do { tmp_o = _grnl(st_BLH_L.omega, o_g2); } while (tmp_o > 4 * o_g2);
			do { tmp_k = _grnl(st_BLH_L.kappa, k_g2); } while (tmp_k > 4 * k_g2);*/
	
		//do { tmp_B = _grnl(st_BLH_R.B, B_g2); } while (tmp_B > 4 * B_g2);
		//do { tmp_L = _grnl(st_BLH_R.L, L_g2); } while (tmp_L > 4 * L_g2);
		//do { tmp_H = _grnl(st_BLH_R.H, H_g2); } while (tmp_H > 4 * H_g2);
		//do { tmp_p = _grnl(st_BLH_R.phi, p_g2); } while (tmp_p > 4 * p_g2);
		//do { tmp_o = _grnl(st_BLH_R.omega, o_g2); } while (tmp_o > 4 * o_g2);
		//do { tmp_k = _grnl(st_BLH_R.kappa, k_g2); } while (tmp_k > 4 * k_g2);

		//////////////////////////////////////////////////////////////////////////
		//位姿的模拟误差应该放在那里合适？
		tmp_B = _grnl(st_BLH_L.B, B_g2);
		tmp_L = _grnl(st_BLH_L.L, L_g2);
		tmp_H = _grnl(st_BLH_L.H, H_g2);
		tmp_p = _grnl(st_BLH_L.phi, p_g2);
		tmp_o = _grnl(st_BLH_L.omega, o_g2);
		tmp_k = _grnl(st_BLH_L.kappa, k_g2);
		Station3DBLH stl(st_BLH_L.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, st_BLH_L.rotate);

		tmp_B = _grnl(st_BLH_R.B, B_g2);
		tmp_L = _grnl(st_BLH_R.L, L_g2);
		tmp_H = _grnl(st_BLH_R.H, H_g2);
		tmp_p = _grnl(st_BLH_R.phi, p_g2);
		tmp_o = _grnl(st_BLH_R.omega, o_g2);
		tmp_k = _grnl(st_BLH_R.kappa, k_g2);

		Station3DBLH str(st_BLH_R.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, st_BLH_R.rotate);

		StaBLH2XYZ(stl, stl_L, a, e2);
		StaBLH2XYZ(str, str_R, a, e2);

		//double ux = 0.0, uy = 0.0, uz = 0.0;
		for (int i = 0; i < point_num; i++)
		{		

			/*do { tmp_a = _grnl(pl_List[i].alpha, px_g2); } while (tmp_a > 4 * px_g2);
			do { tmp_b = _grnl(pl_List[i].beta, py_g2); } while (tmp_b > 4 * py_g2);*/

			tmp_a = _grnl(pl_List[i].alpha, px_g2); tmp_b = _grnl(pl_List[i].beta, py_g2);
			Point2D pl(pl_List[i].imgID, tmp_a, tmp_b);

			/*do { tmp_a = _grnl(pr_List[i].alpha, px_g2); } while (tmp_a > 4 * px_g2);
			do { tmp_b = _grnl(pr_List[i].beta, py_g2); } while (tmp_b > 4 * py_g2);*/

			tmp_a = _grnl(pr_List[i].alpha, px_g2); tmp_b = _grnl(pr_List[i].beta, py_g2);
			Point2D pr(pr_List[i].imgID, tmp_a, tmp_b);

			//前方交会
			vector<Point2D> pointpair;
			vector<Station3D> stationpair;
			pointpair.push_back(pl); pointpair.push_back(pr);
			stationpair.push_back(stl_L); stationpair.push_back(str_R);

			switch (method) {
			case ForwardCross_Linear: {
										  CsinglePointForwardCross sf(pointpair, stationpair);
										  sf.RunSinglePointForwardCross();
										  sf.GetSpacePoint3D(pt_new);
										  X_List[i * max_count + k] = pt_new.X;
										  Y_List[i * max_count + k] = pt_new.Y;
										  Z_List[i * max_count + k] = pt_new.Z;
			}
				break;
			case ForwardCross_Inter: {
										 CintersinglePointForwardCross sisf(pointpair, stationpair,
											 spt_List[i].X, spt_List[i].Y, spt_List[i].Z, 0.001);//初值给10km
										 sisf.RunintersinglePointForwardCross();
										 sisf.GetSpacePoint3D(pt_new);						
										 X_List[i * max_count + k] = pt_new.X;
										 Y_List[i * max_count + k] = pt_new.Y;
										 Z_List[i * max_count + k] = pt_new.Z;
			}
				break;
			default: break;
			}
			/*ux += (pt_new.X - spt_List[i].X) * (pt_new.X - spt_List[i].X);
			uy += (pt_new.Y - spt_List[i].Y) * (pt_new.Y - spt_List[i].Y);
			uz += (pt_new.Z - spt_List[i].Z) * (pt_new.Z - spt_List[i].Z);*/
		}//pointnum

		/*ux = sqrt(ux / point_num); uy = sqrt(uy / point_num); uz = sqrt(uz / point_num);
		printf("i = %d   ux = %lf  uy = %lf uz = %lf\n", k, ux, uy, uz);*/
		//Sleep(1);
	}//count

	/*double *MeanX = CalMean2D(X_List, point_num, max_count);
	double *MeanY = CalMean2D(Y_List, point_num, max_count);
	double *MeanZ = CalMean2D(Z_List, point_num, max_count);*/
	double *Var0X = CalVar2D(X_List, MeanX, point_num, max_count);
	double *Var0Y = CalVar2D(Y_List, MeanY, point_num, max_count);
	double *Var0Z = CalVar2D(Z_List, MeanZ, point_num, max_count);

	printf("计算SEP：\n");
	
	//计算SEP
	//double *SEP = new double[point_num]; memset(SEP, 0, point_num * sizeof(double));
	double *SEP = CalSep2D(X_List, Y_List, Z_List, MeanX, MeanY, MeanZ, point_num, max_count, 100);
	if (fp != NULL) {
		fprintf(fp, "$$ PointID    Ux     Uy      Uz       SEP\n");
		for (int i = 0; i < point_num; i++)
		{
			printf("Var0X : %lf  ", Var0X[i]); printf("Var0Y :  %lf  ", Var0Y[i]); printf("Var0Z : %lf  ", Var0Z[i]);
			printf("SEP   : %lf \n", SEP[i]);
			fprintf(fp, "%5d    %10.6lf   %10.6lf    %10.6lf    %10.6lf\n", i,
				Var0X[i], Var0Y[i], Var0Z[i], SEP[i]);
		}
	}

	delete[] MeanX;
	delete[] MeanY;
	delete[] MeanZ;
	delete[] Var0X;
	delete[] Var0Y;
	delete[] Var0Z;
	delete[] SEP;
}

/*
//此处b代表y方向,l代表x方向
void MetoCaloSimulate(vector<Point2D> pl_List, vector<Point2D> pr_List, vector<SpacePoint3D> spt_List,
	Station3D cmrL, Station3D cmrR, Station3D LeftSP, Station3D RightSP,
	double px_g2, double py_g2, double B_g2, double L_g2, double H_g2,
	double p_g2, double o_g2, double k_g2, int max_count/ * = 5000* /, FILE *fp / *= NULL* /,
	RandomSym sym / *= Normal_Distribution* /, FcrossMethod method / *= ForwardCross_Linear* /)
{
	int point_num = pl_List.size();
	//point_num = 1;

	//Station3D stl_L, str_R;
	SpacePoint3D pt_new;

	double tmp_a, tmp_b;
	double tmp_B, tmp_L, tmp_H;
	double tmp_p, tmp_o, tmp_k;

	double *X_List = new double[point_num * max_count]; memset(X_List, 0, point_num * max_count * sizeof(double));
	double *Y_List = new double[point_num * max_count]; memset(Y_List, 0, point_num * max_count * sizeof(double));
	double *Z_List = new double[point_num * max_count]; memset(Z_List, 0, point_num * max_count * sizeof(double));

	double *MeanX = new double[point_num]; memset(MeanX, 0, point_num * sizeof(double));
	double *MeanY = new double[point_num]; memset(MeanY, 0, point_num * sizeof(double));
	double *MeanZ = new double[point_num]; memset(MeanZ, 0, point_num * sizeof(double));

	for (int i = 0; i < point_num; i++) {
		MeanX[i] = spt_List[i].X; MeanY[i] = spt_List[i].Y; MeanZ[i] = spt_List[i].Z;
	}

	for (int k = 0; k < max_count; k++)
	{
		//////////////////////////////////////////////////////////////////////////
		//位姿的模拟误差应该放在那里合适？
		if (sym == Normal_Distribution) {
			tmp_B = _grnl(LeftSP.Ys, B_g2);
			tmp_L = _grnl(LeftSP.Xs, L_g2);
			tmp_H = _grnl(LeftSP.Zs, H_g2);
			tmp_p = _grnl(LeftSP.phi, p_g2);
			tmp_o = _grnl(LeftSP.omega, o_g2);
			tmp_k = _grnl(LeftSP.kappa, k_g2);
			cmrL.p = tmp_p; cmrL.o = tmp_o; cmrL.k = tmp_k;
			Station3DBLH stl(LeftSP.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, LeftSP.rotate);

			tmp_B = _grnl(RightSP.X, B_g2);
			tmp_L = _grnl(RightSP.L, L_g2);
			tmp_H = _grnl(RightSP.H, H_g2);
			tmp_p = _grnl(RightSP.phi, p_g2);
			tmp_o = _grnl(RightSP.omega, o_g2);
			tmp_k = _grnl(RightSP.kappa, k_g2);
			cmrR.p = tmp_p; cmrR.o = tmp_o; cmrR.k = tmp_k;
			Station3DBLH str(RightSP.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, RightSP.rotate);
		}
		else if (sym == Uniform_Distribution) {
			tmp_B = _rndl(LeftSP.X - 2 * B_g2, LeftSP.X + 2 * B_g2);
			tmp_L = _rndl(LeftSP.L - 2 * L_g2, LeftSP.L + 2 * L_g2);
			tmp_H = _rndl(LeftSP.H - 2 * H_g2, LeftSP.H + 2 * H_g2);
			tmp_p = _rndl(LeftSP.phi - 2 * p_g2, LeftSP.phi + 2 * p_g2);
			tmp_o = _rndl(LeftSP.omega - 2 * o_g2, LeftSP.omega + 2 * o_g2);
			tmp_k = _rndl(LeftSP.kappa - 2 * k_g2, LeftSP.kappa + 2 * k_g2);
			cmrL.p = tmp_p; cmrL.o = tmp_o; cmrL.k = tmp_k;
			Station3DBLH stl(LeftSP.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, LeftSP.rotate);

			tmp_B = _rndl(RightSP.X - 2 * B_g2, RightSP.X + 2 * B_g2);
			tmp_L = _rndl(RightSP.L - 2 * L_g2, RightSP.X + 2 * L_g2);
			tmp_H = _rndl(RightSP.H - 2 * H_g2, RightSP.H + 2 * H_g2);
			tmp_p = _rndl(RightSP.phi - 2 * p_g2, RightSP.phi + 2 * p_g2);
			tmp_o = _rndl(RightSP.omega - 2 * o_g2, RightSP.omega + 2 * o_g2);
			tmp_k = _rndl(RightSP.kappa - 2 * k_g2, RightSP.kappa + 2 * k_g2);
			cmrR.p = tmp_p; cmrR.o = tmp_o; cmrR.k = tmp_k;
			Station3DBLH str(RightSP.imgID, tmp_B, tmp_L, tmp_H, tmp_o, tmp_p, tmp_k, RightSP.rotate);
		}
		else {
			return;
		}

		//StaBLH2XYZ(stl, stl_L, a, e2);
		//StaBLH2XYZ(str, str_R, a, e2);

		//double ux = 0.0, uy = 0.0, uz = 0.0;
		for (int i = 0; i < point_num; i++)
		{
			if (sym == Normal_Distribution) {
				tmp_a = _grnl(pl_List[i].alpha, px_g2); tmp_b = _grnl(pl_List[i].beta, py_g2);
			}
			else if (sym == Uniform_Distribution) {
				tmp_a = _rndl(pl_List[i].alpha - 2 * px_g2, pl_List[i].alpha + 2 * px_g2);
				tmp_b = _rndl(pl_List[i].beta - 2 * py_g2, pl_List[i].beta + 2 * py_g2);
			}
			else {
				return;
			}

			Point2D pl(pl_List[i].imgID, tmp_a, tmp_b);
			if (sym == Normal_Distribution) {
				tmp_a = _grnl(pr_List[i].alpha, px_g2); tmp_b = _grnl(pr_List[i].beta, py_g2);
			}
			else if (sym == Uniform_Distribution) {
				tmp_a = _rndl(pr_List[i].alpha - 2 * px_g2, pr_List[i].alpha + px_g2);
				tmp_b = _rndl(pr_List[i].beta - 2 * py_g2, pr_List[i].beta + py_g2);
			}
			else {
				break;
			}


			Point2D pr(pr_List[i].imgID, tmp_a, tmp_b);

			//前方交会
			vector<Point2D> pointpair;
			vector<Station3DBLH> stationpair;
			vector<CAMERA> cmrs;
			pointpair.push_back(pl); pointpair.push_back(pr);
			cmrs.push_back(cmrL); cmrs.push_back(cmrR);
			stationpair.push_back(LeftSP); stationpair.push_back(RightSP);

			switch (method) {
			case ForwardCross_Linear: {
				CBundleAdjustment cb;
				cb.SetParameters(pointpair, cmrs, stationpair, 2);
				cb.RunForwardCross();
				cb.GetSpacePoint3D(pt_new);
				X_List[i * max_count + k] = pt_new.X;
				Y_List[i * max_count + k] = pt_new.Y;
				Z_List[i * max_count + k] = pt_new.Z;
			}
									  break;
			case ForwardCross_Inter: {
				//CintersinglePointForwardCross sisf(pointpair, stationpair,
				// spt_List[i].X, spt_List[i].Y, spt_List[i].Z, 0.001);//初值给10km
				//sisf.RunintersinglePointForwardCross();
				//sisf.GetSpacePoint3D(pt_new);
				//X_List[i * max_count + k] = pt_new.X;
				//Y_List[i * max_count + k] = pt_new.Y;
				//Z_List[i * max_count + k] = pt_new.Z;
			}
									 break;
			default: break;
			}
		}//pointnum

		//Sleep(1);
	}//count

	double *Var0X = CalVar2D(X_List, MeanX, point_num, max_count);
	double *Var0Y = CalVar2D(Y_List, MeanY, point_num, max_count);
	double *Var0Z = CalVar2D(Z_List, MeanZ, point_num, max_count);

	printf("计算SEP：\n");

	//计算SEP
	double *SEP = new double[point_num]; memset(SEP, 0, point_num * sizeof(double));
	SEP = CalSep2D(X_List, Y_List, Z_List, MeanX, MeanY, MeanZ, point_num, max_count, 100);
	if (fp != NULL) {
		fprintf(fp, "$$ PointID    Ux     Uy      Uz       SEP\n");
		for (int i = 0; i < point_num; i++)
		{
			printf("Var0X : %lf  ", Var0X[i]); printf("Var0Y :  %lf  ", Var0Y[i]); printf("Var0Z : %lf  ", Var0Z[i]);
			printf("SEP   : %lf \n", SEP[i]);
			fprintf(fp, "%5d    %10.6lf   %10.6lf    %10.6lf    %10.6lf\n", i,
				Var0X[i], Var0Y[i], Var0Z[i], SEP[i]);
		}
	}

	delete[] MeanX;
	delete[] MeanY;
	delete[] MeanZ;
	delete[] Var0X;
	delete[] Var0Y;
	delete[] Var0Z;
	delete[] SEP;
}*/