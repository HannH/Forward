#ifndef PROJECT_H
#define PROJECT_H

#include "DataStruct.h"


Point2D invproject(Station3D st, SpacePoint3D pt3d, bool isPlane = false);

void CalRotate(Station3D st, double *R); //������ת����

void POK2R(const double p, const double o, const double k, double *R);
void OPK2R(const double o, const double p, const double k, double *R);

//����������ϵ���������ϵ����ת���� T1 = A * T, T2 = B * T1, A = R(u)*R(i)*R(s)
//     s  - ������ྭ   
//     i  - ������     
//     u  - γ�ȸ���     
//|      | 0  1  0 |         | cos(s)  sin(s)  0 |        |cos(i)  0 -sin(i) |        | cos(u) sin(u) 0 |
//|   B= | 0  0 -1 |  R(s) = |-sin(s)  cos(s)  0 | R(i) = |     0  1       0 | R(u) = |-sin(u) cos(u) 0 |
//|      |-1  0  0 |         |      0       0  1 |        |sin(i)  0  cos(i) |        |      0      0 1 |
//////////////////////////////////////////////////////////////////////////
void CalCoordinateTransfer(const double s, const double i, const double u, double *T);

void CalOPK(const double s, const double i, const double u, const double *R0, double &o, double &p, double &k);
void CalPOK(const double s, const double i, const double u, const double *R0, double &p, double &o, double &k);

#endif