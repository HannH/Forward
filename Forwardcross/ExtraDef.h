#include <string>
typedef struct tagCAMERA{
	double f, fx, fy;
	double u, v;
	double formatx, formaty;
	double sx; //像片横纵比: fx/fy
	double xs, ys, zs;

	double x0, y0, z0;//传感器重心到卫星质心的距离矢量

	double p, o, k;
	double pixelsize;
	RotateSym rotate;
	double k0, k1, k2, k3, p1, p2, b1, b2;

	double R0; //畸变改正零点距离
	int flag;				//相机主点改正标志，对应于主点改正的正负号。原始检校结果对应0，每逆时针旋转90度加1，即第二象限为1，即第三象限为2，第四象限为3

	int CameraID;
	int imageID; //影像编号,从零开始的唯一序列号
	int m_ImageName;	//数字形式的影像编号，可以任意编号，无需从零开始
	char ImageName[256];  //xxw, 2015-08-29,存储bai文件中的影像名
	double wtf, wtu, wtv;
	double wtxs, wtys, wtzs;
	double wtphi, wtomega, wtkappa;
	tagCAMERA()
	{
		strcpy(ImageName, "");      //2015-08-29,xxw

		f = fx = fy = 1;
		u = v = 0;
		formatx = formaty = 0;
		sx = 1;
		xs = ys = zs = 0;
		x0 = y0 = z0 = 0;
		p = o = k = 0;
		rotate = RotateSym::pok;
		pixelsize = 0.008;
		k0 = k1 = k2 = k3 = p1 = p2 = b1 = b2 = 0;
		R0 = 0;
		flag = 0;
		CameraID = 0;
		imageID = 0;

		m_ImageName = 0;

		wtf = wtu = wtv = 0;
		wtxs = wtys = wtzs = 0;
		wtphi = wtomega = wtkappa = 0;
	}
}CAMERA;
