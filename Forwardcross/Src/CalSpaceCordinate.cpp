#include <stdio.h>
#include <stdlib.h>

#include "DataSimulater.h"
#include "SinglePointFcoss.h"
#include "SGCordinate.h"
#include "Project.h"
#include "Matrix.h"

#pragma warning (disable : 4018)
#pragma warning (disable : 4101)

/*
int main(int argc, char **argv)
{
	double B1 = 0, L1 = 100, H1 = 35786;//高度36000km
	double B2 = 0, L2 = 130, H2 = 35786;
	Cg cg(6378.137, 0.0066943799013);
	
	double Xs1, Ys1, Zs1, Xs2, Ys2, Zs2;
	cg.g2r(B1, L1, H1, Xs1, Ys1, Zs1);
	cg.g2r(B2, L2, H2, Xs2, Ys2, Zs2);

	printf("ST1: %lf %lf %lf\nST2: %lf %lf %lf\n", Xs1, Ys1, Zs1, Xs2, Ys2, Zs2);

	double p1, o1, k1, p2, o2, k2;
	double R1[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }, R2[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

	CalPOK(-L1 * PI / 180, -0 * PI / 180, 0,  R1, p1, o1, k1);
	CalPOK(-L2 * PI / 180, -0 * PI / 180, 0,  R2, p2, o2, k2);

	Station3D st1(0, Xs1, Ys1, Zs1, o1, p1, k1, pok); //姿态未统一坐标系
	Station3D st2(1, Xs2, Ys2, Zs2, o2, p2, k2, pok);

	double B3 = 5, L3 = 115, H3 = 500; //高度500km
	double X, Y, Z;
	cg.g2r(B3, L3, H3, X, Y, Z);

	/ *X = -2533.949447;
	Y = 5434.073204;
	Z = 0.207178;* /
	
	printf("目标坐标： %lf %lf %lf\n", X, Y, Z);

	SpacePoint3D spt(X, Y, Z);
	
	Point2D ptl = invproject(st1, spt);
	Point2D ptr = invproject(st2, spt);

	printf("alpha = %.20lf beta = %.20lf\nalpha = %.20lf beta = %.20lf\n",
		ptl.alpha , ptl.beta, ptr.alpha, ptr.beta);

	Point2D pt_1(0, ptl.alpha - 60e-6, ptl.beta - 60e-6);
	Point2D pt_2(1, ptr.alpha - 120e-6, ptr.beta - 30e-6);

	vector<Point2D> vpt;
	vector<Station3D> vst;

	vpt.push_back(pt_1); vpt.push_back(pt_2);
	vst.push_back(st1);	vst.push_back(st2);

	printf("目标坐标： %lf %lf %lf\n", spt.X, spt.Y, spt.Z);
	//CintersinglePointForwardCross sitspf(vpt, vst, spt.X+10, spt.Y-10, spt.Z+2, 0.01);
	//sitspf.RunintersinglePointForwardCross();
	//double m0, Qxx[9];
	//sitspf.GetSpacePoint3D_m0(m0);
	//sitspf.GetSpacePoint3D_Qxx(Qxx);

	//printf("单位权方差： %.20lf\n", m0);
	//printf("精度： Dxx = %lf Dyy = %lf Dzz = %lf\n", / *m0 ** / Qxx[0], / *m0 * * /Qxx[4], / *m0 * * /Qxx[8]);

	CsinglePointForwardCross spf(vpt, vst);
	spf.RunSinglePointForwardCross();

	SpacePoint3D pt3d;
	spf.GetSpacePoint3D(pt3d);

	printf("前方交会(线性) ： X = %lf\n Y = %lf\n Z = %lf\n", pt3d.X, pt3d.Y, pt3d.Z);

	/ *sitspf.GetSpacePoint3D(pt3d);
	printf("前方交会（迭代） ： X = %lf\n Y = %lf\n Z = %lf\n", pt3d.X, pt3d.Y, pt3d.Z);* /

	double B, L, H;
	cg.r2g(pt3d.X, pt3d.Y, pt3d.Z, B, L, H);
	printf("经纬高 ： B = %lf\n L = %lf\n H = %lf\n", B, L, H);

	cg.g2r(B, L, H, spt.X, spt.Y, spt.Z);
	ptl = invproject(st1, spt);
	ptr = invproject(st2, spt);

	printf("alpha = %lf beta = %lf\nalpha = %lf beta = %lf",
		ptl.alpha * RAD2DEG, ptl.beta* RAD2DEG, ptr.alpha* RAD2DEG, ptr.beta* RAD2DEG);

	system("pause");
	return 0;

}*/

int main1(int argc, char **argv)
{
	char *xnfile = "C:\\Users\\xxw\\Desktop\\0SimulationPoints.txt";
	vector<WxVirtualPoint> VirtualPoint;
	LoadVirtualPoint(xnfile, VirtualPoint);

	Point2D pl = VirtualPoint[0].pl; pl.imgID = 0;
	Point2D pr = VirtualPoint[0].pr; pr.imgID = 1;

	printf("%lf %lf %lf %lf\n", pl.alpha, pl.beta, pr.alpha, pr.beta);

	//产生模拟数据
	double threshold_px = 120e-6;
	int threshold_num = 5000/*10*//*10000*/;
	vector<Point2D> vpl_List, vpr_List;
	
	printf("开始模拟数据！\n");
	int sim_num = SimulatePointPairList(pl, pr, vpl_List, vpr_List, threshold_px, threshold_num);
	if ((double)sim_num / threshold_num < 0.5) {
		printf("模拟的点数太少！  %d\n", sim_num);
		system("pause");
		exit(0);
	}

	printf("模拟数目：n1 = %d, n2 = %d\n", vpl_List.size(), vpr_List.size());

	//for (int i = 0; i < sim_num; i++)
	//{
	//	printf("lid = %d, rid = %d\n", vpl_List[i].imgID, vpr_List[i].imgID);
	//}
	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);

	Station3D Sta_L, Sta_R;
	StaBLH2XYZ(StaBLH_L, Sta_L, WGS84_a, WGS84_e2);
	StaBLH2XYZ(StaBLH_R, Sta_R, WGS84_a, WGS84_e2);

	printf("开始前方交会!\n");
	vector<SpacePoint3D> vpt3d;
	MulPointSingleForwardCross(vpl_List, vpr_List, Sta_L, Sta_R, vpt3d);

	printf("开始保存结果!\n");
	char *savefile = "C:\\Users\\xxw\\Desktop\\result.txt";
	SaveSpacePoint3D(savefile, vpt3d);

	ApproximateSEP(vpt3d, 5);

	ApproximateSEP2(vpt3d, VirtualPoint[0].pt3d, 5);

	system("pause");
	return 0;
}

int main2(int argc, char **argv)
{
	const char *xnfile = "C:\\Users\\xxw\\Desktop\\0SimulationPoints.txt";
	vector<WxVirtualPoint> VirtualPoint;
	LoadVirtualPoint(xnfile, VirtualPoint);

	const char *sepfile = "C:\\Users\\xxw\\Desktop\\SepResult正太.txt";

	FILE *fp = fopen(sepfile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}
	fprintf(fp, "$$ SEP Number\n");
	fprintf(fp, "$$ ID  SEP\n");
	fprintf(fp, "%d\n", sepfile);

	double threshold_px = 120e-6;
	int threshold_maxnum = 5000/*10*//*10000*/;
	int threshold_minnum = 50;
	int threshold_calNum = 1000;
	double sep;

	Point2D pl, pr;
	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);

	for (int i = 0; i < /*1*/VirtualPoint.size(); i++)
	{
		printf("开始算第 %d 个点的SEP！\n", i + 1);
		pl = VirtualPoint[i].pl; pl.imgID = 0;
		pr = VirtualPoint[i].pr; pr.imgID = 1;

		sep = PointPairLocationSep(pl, pr, VirtualPoint[i].pt3d, StaBLH_L, StaBLH_R, threshold_px, threshold_maxnum); //单次模拟
		/*sep = StatisticPointPairSEP(pl, pr, VirtualPoint[i].pt3d, StaBLH_L, StaBLH_R,
			threshold_px, threshold_maxnum, threshold_minnum, threshold_calNum, WGS84_a, WGS84_e2, Normal_Distribution);*/

		fprintf(fp, "%4d   %10.6lf\n", i, sep);
		printf("sep[%d] : %lf\n", i, sep);
	}
	fclose(fp);

	system("pause");

	return 0;
}

int main3(int argc, char **argv)
{
	const char *xnfile = "C:\\Users\\xxw\\Desktop\\0SimulationPoints.txt";
	vector<WxVirtualPoint> VirtualPoint;
	LoadVirtualPoint(xnfile, VirtualPoint);

	const char *sepfile = "C:\\Users\\xxw\\Desktop\\InterMaxResult.txt";

	FILE *fp = fopen(sepfile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}
	fprintf(fp, "$$ MaxResidual Number\n");
	fprintf(fp, "$$ ID  dX   dY   dZ  REMS\n");
	fprintf(fp, "%d\n", VirtualPoint.size());

	double max_threshold_px = 120e-6;
	double sep;

	Point2D pl, pr;
	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);
	Station3D Sta_L, Sta_R;
	StaBLH2XYZ(StaBLH_L, Sta_L, WGS84_a, WGS84_e2);
	StaBLH2XYZ(StaBLH_R, Sta_R, WGS84_a, WGS84_e2);

	double dX[3] = { 0.0 };

	for (int i = 0; i < /*1*/VirtualPoint.size(); i++)
	{
		printf("开始算第 %d 个点的最大残差！\n", i + 1);
		pl = VirtualPoint[i].pl; pl.imgID = 0;
		pr = VirtualPoint[i].pr; pr.imgID = 1;

		sep = MaxResidulPointPair(pl, pr, Sta_L, Sta_R, VirtualPoint[i].pt3d, max_threshold_px, dX, ForwardCross_Inter);
		
		fprintf(fp, "%4d   %10.6lf  %10.6lf  %10.6lf   %10.6lf\n", i, dX[0], dX[1], dX[2], sep);
		printf("sep[%d] : %lf\n", i, sep);
	}
	fclose(fp);

	system("pause");

	return 0;
}

int main4(int argc, char **argv)
{
	const char *xnfile = "C:\\Users\\xxw\\Desktop\\0SimulationPoints.txt";
	vector<WxVirtualPoint> VirtualPoint;
	LoadVirtualPoint(xnfile, VirtualPoint);

	const char *sepfile = "C:\\Users\\xxw\\Desktop\\InterMaxResult蒙特卡洛.txt";

	FILE *fp = fopen(sepfile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}
	fprintf(fp, "$$ MaxResidual Number\n");
	fprintf(fp, "$$ ID  dX   dY   dZ  REMS\n");
	fprintf(fp, "%d\n", VirtualPoint.size());

	double max_threshold_px = 120e-6;
	double max_B = 0e-6;
	double max_L = 0e-6;
	double max_H = 0e-3;
	double max_p = 3.4 / 3600/*20e-6*/;//资三1.7弧秒的精度
	double max_o = 3.4 / 3600/*20e-6*/;
	double max_k = 3.4 / 3600/*20e-6*/;
	int max_num = 10000;
	double sep;

	Point2D pl, pr;
	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);
	Station3D Sta_L, Sta_R;
	StaBLH2XYZ(StaBLH_L, Sta_L, WGS84_a, WGS84_e2);
	StaBLH2XYZ(StaBLH_R, Sta_R, WGS84_a, WGS84_e2);

	//double dX[3] = { 0.0 };
	double *dX = new double[3];
	memset(dX, 0, 3 * sizeof(double));

	for (int i = 0; i < /*1*/VirtualPoint.size(); i++)
	{
		printf("开始算第 %d 个点的残差！\n", i + 1);
		pl = VirtualPoint[i].pl; pl.imgID = 0;
		pr = VirtualPoint[i].pr; pr.imgID = 1;

		dX = SymstemResidualForwardCross(pl, pr, StaBLH_L, StaBLH_R, VirtualPoint[i].pt3d,
			max_threshold_px, max_threshold_px, max_B, max_H, max_p, max_o, max_k, max_num);

		sep = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

		fprintf(fp, "%4d   %10.6lf  %10.6lf  %10.6lf   %10.6lf\n", i, dX[0], dX[1], dX[2], sep);
		printf("sep[%d] : %lf\n", i, sep);
	}
	fclose(fp);
	delete[] dX;

	system("pause");

	return 0;
}


//基本确定，改了转角公式
int main5(int argc, char **argv)
{
	const char *xnfile = "C:\\Users\\xxw\\Desktop\\0SimulationPoints.txt";
	vector<SpacePoint3D> spt_List;
	vector<Point2D> pl_List, pr_List;

	printf("读取虚拟点文件！\n");
	LoadVirtulPoint2File(xnfile, spt_List, pl_List, pr_List);
	printf("共有点数：%d\n", spt_List.size());

	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);

	double px_g2 = 30e-6, py_g2 = 30e-6; //4倍中误差 = 120urad
	double B_g2 = 100e-6, L_g2 = 100e-6, H_g2 = 30e-3;
	double p_g2 = 50e-6 /*1.7/ 3600*/, o_g2 = 50e-6 /*1.7/ 3600*/, k_g2 = 50e-6/*1.7/ 3600*/;
	int max_count = 10000;

	const char *sepfile = "C:\\Users\\xxw\\Desktop\\LinerSEP蒙特卡洛.txt";
	FILE *fp = fopen(sepfile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}

	printf("开始MetoCalo模拟!\n");
	MetoCaloSimulate(pl_List, pr_List, StaBLH_L, StaBLH_R, spt_List, px_g2, py_g2,
		B_g2, L_g2, H_g2, p_g2, o_g2, k_g2, max_count, fp);

	fclose(fp);
	system("pause");
	return 0;
}

int main6(int argc, char **argv)
{
	Station3DBLH StaBLH_L(0, 0, 100, 36000, 0, 0, 0, pok); //标称图的位姿
	Station3DBLH StaBLH_R(1, 0, 130, 36000, 0, 0, 0, pok);

	Station3D Sta_L, Sta_R;
	StaBLH2XYZ(StaBLH_L, Sta_L, WGS84_a, WGS84_e2);
	StaBLH2XYZ(StaBLH_R, Sta_R, WGS84_a, WGS84_e2);

	double R[9], R0[9];
	POK2R(Sta_L.phi, Sta_L.omega, Sta_L.kappa, R);
	printf("R:\n");
	CMatrix matrixapp;
	matrixapp.MatrixPrint(R, 3, 3);

	R0[0] =  -sin(-100 * DEG2RAD); R0[1] = 0; R0[2] = -cos(-100 * DEG2RAD);
	R0[3] =  cos(-100 * DEG2RAD); R0[4] = 0; R0[5] = sin(-100 * DEG2RAD);
	R0[6] = 0; R0[7] = 1; R0[8] = 0;

	matrixapp.MatrixPrint(R0, 3, 3);
 
	system("pause");
	return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>

#include "..\ExtraDef.h"
#include "..\Src\DataStruct.h"
#include "..\Src\SGCordinate.h"
#include "..\Src\Project.h"
#include "..\Src\Matrix.h"
#include "DataSimulater.h"
#include "..\Fcross.h"
const double pi=3.14159265357 ;
#pragma warning (disable : 4018)
#pragma warning (disable : 4101)

int main(int argc, char **argv)
{
	//Station3DBLH stl(0, 0, 000, 35786, 0, 0, 0, pok);
	//Station3DBLH str(1, 0, 100, 35786, 0, 0, 0, pok);
	Station3D LeftCamera(3000.0,4000.0,2886.7513,pi/3.0,0,0,pok)
		, RightCamera(3402.9847, 3663.2902, 2886.7513, pi / 3.0, 0, 0, pok);



	double X=100, Y=100, Z=50;

	SpacePoint3D spt(X, Y, Z);

	Point2D ptl = invproject(LeftCamera,spt,true);
	Point2D ptr = invproject(RightCamera,spt,true);

	//Point2D pt2 = Point2AnotherImg(ptl, stl, str, cmrl, cmrr);

	//printf("坐标转化后： alpha = %lf,  beta = %lf\n", pt2.alpha*RAD2DEG, pt2.beta*RAD2DEG);

	printf("lx = %.20lf ly = %.20lf\nrx = %.20lf ry = %.20lf\n",
		ptl.x , ptl.y , ptr.x , ptr.y );

	//前方交会
	Point2D pt_1(0, ptl.alpha /*- 60e-6*/, ptl.beta/* - 60e-6*/);
	Point2D pt_2(1, ptr.alpha /*- 120e-6*/, ptr.beta /*- 30e-6*/);
	vector<Point2D> vpt2d; vpt2d.push_back(pt_1); vpt2d.push_back(pt_2);

	vector<Station3D>  CameraVec;  CameraVec.push_back(LeftCamera); CameraVec.push_back(RightCamera);
	vector<SpacePoint3D> StationVec; StationVec.push_back(spt);

	CBundleAdjustment cb;
	cb.SetInnerEle(0, 0, -0.1);
	cb.SetParameters(vpt2d, CameraVec, StationVec, 2);
	cb.RunForwardCross();

	SpacePoint3D pt3d;
	cb.GetSpacePoint3D(pt3d);

	printf("坐标： X = %lf, Y = %lf, Z = %lf\n", pt3d.X, pt3d.Y, pt3d.Z);
	//掺正态误差
	double delta = 0.2;
	pt_1.x=_grnl(pt_1.x, delta);
	pt_1.y = _grnl(pt_1.y, delta);
	pt_2.x = _grnl(pt_2.x, delta);
	pt_2.y = _grnl(pt_2.y, delta);


	/*//SEP
	vector<Point2D> PL_List, PR_List;
	PL_List.push_back(ptl); PR_List.push_back(ptr);

	double px_g2 = 60e-6, py_g2 = 60e-6; //2倍中误差 = 120urad
	double B_g2 = 100e-6, L_g2 = 100e-6, H_g2 = 30e-3;
	double p_g2 = 50e-6 / *0.0/ 3600* /, o_g2 = 50e-6 / *0.0/ 3600* /, k_g2 = 50e-6/ *0.0/ 3600* /;
	int max_count = 10000;

	vector<SpacePoint3D> vpt3d; vpt3d.push_back(spt);

	const char *sepfile = "C:\\Users\\xxw\\Desktop\\New_LinerSEP蒙特卡洛.txt";
	FILE *fp = fopen(sepfile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}

	MetoCaloSimulate(PL_List, PR_List, vpt3d, cmrl, cmrr, stl, str,
		px_g2, py_g2, B_g2, L_g2, H_g2, p_g2, o_g2, k_g2, max_count,
		fp, Normal_Distribution, ForwardCross_Linear);
	fclose(fp);*/

	system("pause");
	return 0;
}

/*
int main(int argc, char **argv)
{
	char *GDFile;
	char SEPFile[256];
	GDFile = / *"C:\\Users\\xxw\\Desktop\\数据模拟\\Li\\相向\\Simulation_points_NoDT_0.txt"* /argv[1];

	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];

	_splitpath_s(GDFile, drive, dir, fname, ext);
	sprintf_s(SEPFile, "%s%s%s_SEP_位置.txt", drive, dir, fname);

	vector<Point2D> Lpt2d_List, Rpt2d_List;
	vector<SpacePoint3D> Lspt3d_List, Rspt3d_List;

	LoadVirtulPoint2File(GDFile, Lspt3d_List, Rspt3d_List, Lpt2d_List, Rpt2d_List);

	Station3DBLH stl(0, 0, 100, 35786, 0, 0, 0, pok);
	Station3DBLH str(1, 0, 140, 35786, 0, 0, 0, pok);

	CAMERA cmrl, cmrr;
	cmrl.u = cmrl.v = cmrr.u = cmrr.v = 0;
	cmrr.p = cmrr.o = cmrr.k = cmrl.p = cmrl.o = cmrl.k = 0.0;

	double px_g2 = 60e-6, py_g2 = 60e-6; //2倍中误差 = 120urad
	double B_g2 = 0/ *0.5 / 3600* / / * / 3600* // *000e-6* /, L_g2 = 0 / * / 3600* // *000e-6* /, H_g2 = 000e-3;
	double p_g2 = / *100e-6* / 0.5 / 3600, o_g2 = / *100e-6* / 0.5 / 3600, k_g2 = / *100e-6* /0.5 / 3600;
	int max_count = 10000;

	FILE *fp = fopen(SEPFile, "w");
	if (!fp)
	{
		printf("SEP文件保存失败!\n");
		return 0;
	}

	MetoCaloSimulate(Lpt2d_List, Rpt2d_List, Lspt3d_List, cmrl, cmrr, stl, str,
		px_g2, py_g2, B_g2, L_g2, H_g2, p_g2, o_g2, k_g2, max_count,
		fp, Normal_Distribution, ForwardCross_Linear, WGS84_a, WGS84_e2);

	fclose(fp);

	return 0;
}*/