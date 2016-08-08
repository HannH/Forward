#include <string>
typedef struct tagCAMERA{
	double f, fx, fy;
	double u, v;
	double formatx, formaty;
	double sx; //��Ƭ���ݱ�: fx/fy
	double xs, ys, zs;

	double x0, y0, z0;//���������ĵ��������ĵľ���ʸ��

	double p, o, k;
	double pixelsize;
	RotateSym rotate;
	double k0, k1, k2, k3, p1, p2, b1, b2;

	double R0; //�������������
	int flag;				//������������־����Ӧ����������������š�ԭʼ��У�����Ӧ0��ÿ��ʱ����ת90�ȼ�1�����ڶ�����Ϊ1������������Ϊ2����������Ϊ3

	int CameraID;
	int imageID; //Ӱ����,���㿪ʼ��Ψһ���к�
	int m_ImageName;	//������ʽ��Ӱ���ţ����������ţ�������㿪ʼ
	char ImageName[256];  //xxw, 2015-08-29,�洢bai�ļ��е�Ӱ����
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
