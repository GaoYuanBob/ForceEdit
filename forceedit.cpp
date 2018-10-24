/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
History
$Log: meshedit.cpp,v $
****************************************************************************/
#include <QtGui>
#include <stdio.h>
#include <vector>

#include "common/curves.h"
#include "common/realtimebuildgum.h"
#include "forceedit.h"
#include <vcg/complex/trimesh/append.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <gl/glut.h>
#include <Eigen/Dense>
#include <wrap/gl/pick.h>
#include <wrap/qt/gl_label.h>
#include <wrap/gui/trackball.h>
#include <wrap/gui/coordinateframe.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_stl.h>
#include <wrap/io_trimesh/import_stl.h>

///////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////// - GY
#define Kr 1.0000  //牙齿的刚度系数
#define Kc 10.000 //牙套的刚度系数

#include <QtGui>

#include <math.h>
#include <CString>
#include <cstring>
#include <stdlib.h>
#include <meshlab/glarea.h>
#include "forceedit.h"
#include <wrap/gl/pick.h>
#include <vcg/complex/trimesh/append.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/qt/gl_label.h>
#include "common/Attachmodel.h"
#include <GL/glut.h>
#include <direct.h>//C++创建文件夹相关头文件
#include <fstream> 
#include <sstream>  
#include <Windows.h>
#include <QFileInfo>
#include <Eigen/Dense>
#include <vcg/complex/trimesh/create/platonic.h>
#include <QProgressBar>

#include "CollisionSimulate.h"
#include "CSG/files.h"
#include "GetUndercut.h"

#include <QFileInfoList>  
#include <QDir>
#include <string>
#include <windows.h>
#include <string>
#include <fstream>

// 用来牙套裁剪
#include "GeodeticCalculator_Dijk.h"
#include "DijkstraSet.h"
#include "Mesh.h"
#include "common/interfaces.h"
/////////////////////////////////////////////////////// GY
using namespace std;
using namespace vcg;
using namespace Eigen;

const float FPI = 3.141592653;
const float FPi = 3.141592653;
#define Kr 1.0000
#define L 3 
#define lamda 0
std::vector<ToothNode* > RTBGum::newseq;
std::vector<ToothNode* > RTBGum::meshTreeAll;
std::vector<Pontic* > RTBGum::pontics;
vcg::Point3f RTBGum::planeN;

std::vector<NearestPoint> RTBGum::PointsPairs;
std::vector<vcg::Point3f> RTBGum::gumPoints;//按一定规则存放所有牙龈线控制点
std::vector<vcg::Point3f> RTBGum::gumPointsBottom;//存放投影到底部的边界线控制点
std::vector<vcg::Point3f> RTBGum::gumPointsBottomPart;//存放投影到底部控制点
std::vector<vcg::Point3f> RTBGum::gumFacePoints;//存放牙龈侧面的顶点
std::vector<vcg::Point3f> RTBGum::gumFacePointsTop;//存放牙龈两端的顶点（头）
std::vector<vcg::Point3f> RTBGum::gumFacePointsTail;//存放牙龈两端的顶点（尾）
std::vector<vcg::Point3f> RTBGum::gumFacePointsBL;//存放牙龈侧面的顶点(颊舌侧)

std::vector<vcg::Point3f> RTBGum::safcPoints;//存放底部拟合外侧曲线
std::vector<vcg::Point3f> RTBGum::safcPoints2;//存放底部拟合内侧曲线
std::vector<vcg::Point3f> RTBGum::controlPoints;//存放底部外侧控制点
std::vector<vcg::Point3f> RTBGum::controlPointsDir;//存放底部外侧控制点方向
std::vector<vcg::Point3f> RTBGum::controlPoints2;//存放底部内侧控制点
std::vector<vcg::Point3f> RTBGum::controlPointsDir2;//存放底部内侧控制点方向
std::vector<vcg::Point3f> RTBGum::controlPointsOrigin;//存放底部外侧控制点初始值
std::vector<vcg::Point3f> RTBGum::controlPointsOrigin2;//存放底部内侧控制点初始值

std::vector<vcg::Point3f> RTBGum::controlPointsHeadTopTemp;//存放两端的控制点（头）
std::vector<vcg::Point3f> RTBGum::controlPointsTailTopTemp;//存放两端的控制点（尾）


std::vector<vcg::Point3f> RTBGum::gumFacePointsTailTemp;//底部曲线头
std::vector<vcg::Point3f> RTBGum::gumFacePointsTailTemp2;//底部曲线尾

std::vector<vcg::Point3f> RTBGum::teethBordersInOrderBottom;//存放投影到底部的有序边界点

static std::vector<vcg::Point3f> CSs;

ForceEditPlugin::ForceEditPlugin() : rubberband(Color4b(255,170,85,255))
{
	forceJudge=false;
	forceui = 0;
	has_init = false;
	pickloc = false;
	isRightClicked = false;
	hasnewfata = false;
	isMoving = false;
	was_ready = false;
	is_measure = false;
	firstclick = false;
	current_teeth = 0;
	curMode = ManuMode::None;
	attype = 0;
	importID = -1;
	inputValue = 0.0;
	original_Transform = vcg::Matrix44f::Identity();
	delta_Transform = vcg::Matrix44f::Identity();
}

const QString ForceEditPlugin::Info() 
{
	return tr("Force Edit.");
}

inline Point3f M(Matrix44f res,Point3f v)
{
	vcg::Point3f Vn;
	Vn.X()=res[0][0]*v.X()+res[0][1]*v.Y()+res[0][2]*v.Z();
	Vn.Y()=res[1][0]*v.X()+res[1][1]*v.Y()+res[1][2]*v.Z();
	Vn.Z()=res[2][0]*v.X()+res[2][1]*v.Y()+res[2][2]*v.Z();
	return Vn;
}

inline Point3f M(float res[3][3],Point3f v)
{
	vcg::Point3f Vn;
	Vn.X()=res[0][0]*v.X()+res[0][1]*v.Y()+res[0][2]*v.Z();
	Vn.Y()=res[1][0]*v.X()+res[1][1]*v.Y()+res[1][2]*v.Z();
	Vn.Z()=res[2][0]*v.X()+res[2][1]*v.Y()+res[2][2]*v.Z();
	return Vn;
}

inline Point3f M(double res[3][3],Point3f v)
{
	vcg::Point3f Vn;
	Vn.X()=res[0][0]*v.X()+res[0][1]*v.Y()+res[0][2]*v.Z();
	Vn.Y()=res[1][0]*v.X()+res[1][1]*v.Y()+res[1][2]*v.Z();
	Vn.Z()=res[2][0]*v.X()+res[2][1]*v.Y()+res[2][2]*v.Z();
	return Vn;
}

inline Point3f M(Point3f ceshipoints,Point3f m,Point3f x,Point3f y,Point3f z){//¼ÆËãÏòÁ¿v¾­¹ý¾ØÕóMµÄ±ä»»ºóµÄÖµ
	vcg::Point3f Vn;
	Vn.X()=ceshipoints.X()*x.X()+ceshipoints.Y()*x.Y()+ceshipoints.Z()*x.Z()-(m.X()*x.X()+m.Y()*x.Y()+m.Z()*x.Z());
	Vn.Y()=ceshipoints.X()*y.X()+ceshipoints.Y()*y.Y()+ceshipoints.Z()*y.Z()-(m.X()*y.X()+m.Y()*y.Y()+m.Z()*y.Z());
	Vn.Z()=ceshipoints.X()*z.X()+ceshipoints.Y()*z.Y()+ceshipoints.Z()*z.Z()-(m.X()*z.X()+m.Y()*z.Y()+m.Z()*z.Z());
	return Vn;
}

//---------------------------------------------------------------------------------------
void ForceEditPlugin::DrawWireCone(vcg::Point3f loc,vcg::Point3f direct)
{
	vcg::Point3f di = direct.Normalize();
	vcg::Point3f zasix(0,0,1);
	vcg::Point3f roasix = zasix ^ di;
	roasix.Normalize();
	double rotavalue = acos((di * zasix)/direct.Norm())/FPi*180;
	glColor4f(1.0f,1.0f,0.0f,0.3f);

	glPushMatrix();
	glTranslatef(loc.X(), loc.Y(), loc.Z());
	glRotatef(rotavalue,roasix[0],roasix[1],roasix[2]);
	glPushMatrix();
	glTranslatef(0, 0, -10);
	//画线框圆锥
	glutWireCone(4,10,32,1);
	glColor4f(0.0f,0.0f,1.0f,0.3f);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//画底面原盘
	GLUquadricObj *objCylinder1 = gluNewQuadric();
	gluDisk(objCylinder1,0,4,32,5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, 0, -5);
	//画中间原盘
	GLUquadricObj *objCylinder2 = gluNewQuadric();
	gluDisk(objCylinder2,0,2,32,5);
	glPopMatrix();
	glPopMatrix();
	glDisable(GL_BLEND);

}

void ForceEditPlugin::DrawVector(vcg::Point3f &vector,vcg::Color4b color)
{
	float beishu = 1;
	vcg::Point3f zasix(0,0,1);
	vcg::Point3f roasix = zasix ^ vector;
	roasix.Normalize();
	double rotavalue = acos((vector * zasix)/vector.Norm())/FPi*180;
	vcg::Point3f center_p = CMN->barycentric_coord;

	glPushMatrix();
	glTranslatef(center_p.X(), center_p.Y(), center_p.Z());
	glRotatef(rotavalue,roasix[0],roasix[1],roasix[2]);
	glColor3f(color[0]/255.0, color[1]/255.0, color[2]/255.0);

	GLUquadricObj *objCylinder3 = gluNewQuadric();
	gluDisk(objCylinder3,0,0.5,32,5);

	GLUquadricObj *objCylinder1 = gluNewQuadric();
	gluCylinder(objCylinder1, 0.3, 0.3, beishu * vector.Norm(), 32, 5);

	glPushMatrix();
	glTranslatef(0,0,beishu * vector.Norm());
	glutSolidCone(0.5,2,32,5);
	GLUquadricObj *objCylinder2 = gluNewQuadric();
	gluDisk(objCylinder2,0,0.3,32,5);
	glPopMatrix();
	glPopMatrix();
}

//---------------------------------------------------------------------------------------
void ForceEditPlugin::DrawCubes(float r, float g, float b)
{
	glColor4f(r,g,b,1.0);
	glBegin (GL_LINES);
	// mid line
	glVertex3f( 0.0,  0.0, -1.0);
	glVertex3f( 0.0,  0.0,  1.0);
	glEnd ();

	glBegin (GL_LINES);
	// right cube
	glVertex3f(  0.0,  0.0,  1.0);
	glVertex3f(  0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.0,  1.0);
	glVertex3f( -0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.0,  1.0);
	glVertex3f(  0.0, -0.1,  1.1);
	glVertex3f(  0.0,  0.0,  1.0);
	glVertex3f(  0.0,  0.1,  1.1);
	glVertex3f(  0.0,  0.0,  1.2);
	glVertex3f(  0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.0,  1.2);
	glVertex3f( -0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.0,  1.2);
	glVertex3f(  0.0, -0.1,  1.1);
	glVertex3f(  0.0,  0.0,  1.2);
	glVertex3f(  0.0,  0.1,  1.1);
	glEnd ();

	glBegin (GL_LINES);
	// left cube
	glVertex3f(  0.0,  0.0, -1.0);
	glVertex3f(  0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.0, -1.0);
	glVertex3f( -0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.0, -1.0);
	glVertex3f(  0.0, -0.1, -1.1);
	glVertex3f(  0.0,  0.0, -1.0);
	glVertex3f(  0.0,  0.1, -1.1);
	glVertex3f(  0.0,  0.0, -1.2);
	glVertex3f(  0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.0, -1.2);
	glVertex3f( -0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.0, -1.2);
	glVertex3f(  0.0, -0.1, -1.1);
	glVertex3f(  0.0,  0.0, -1.2);
	glVertex3f(  0.0,  0.1, -1.1);
	glEnd ();


	// right cube
	glColor4f(std::min(1.0f,r+0.2f), std::min(1.0f,g+0.2f), std::min(1.0f,b+0.2f),0.5);
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0,  1.2);
	glVertex3f(  0.0,  0.1,  1.1);
	glVertex3f( -0.1,  0.0,  1.1);
	glVertex3f(  0.0, -0.1,  1.1);
	glVertex3f(  0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.1,  1.1);
	glEnd();
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0,  1.0);
	glVertex3f(  0.0,  0.1,  1.1);
	glVertex3f( -0.1,  0.0,  1.1);
	glVertex3f(  0.0, -0.1,  1.1);
	glVertex3f(  0.1,  0.0,  1.1);
	glVertex3f(  0.0,  0.1,  1.1);
	glEnd();

	// left cube
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0, -1.2);
	glVertex3f(  0.0,  0.1, -1.1);
	glVertex3f( -0.1,  0.0, -1.1);
	glVertex3f(  0.0, -0.1, -1.1);
	glVertex3f(  0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.1, -1.1);
	glEnd();
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0, -1.0);
	glVertex3f(  0.0,  0.1, -1.1);
	glVertex3f( -0.1,  0.0, -1.1);
	glVertex3f(  0.0, -0.1, -1.1);
	glVertex3f(  0.1,  0.0, -1.1);
	glVertex3f(  0.0,  0.1, -1.1);
	glEnd();
}

void ForceEditPlugin::DrawArrows(float r, float g, float b)
{
	glColor4f(r,g,b,1.0);
	glBegin (GL_LINES);
	// mid line
	glVertex3f( 0.0,  0.0, -1.1);
	glVertex3f( 0.0,  0.0,  1.1);

	// right arrow
	glVertex3f(  0.0,  0.0,  1.1);
	glVertex3f(  0.1,  0.1,  0.9);
	glVertex3f(  0.0,  0.0,  1.1);
	glVertex3f( -0.1,  0.1,  0.9);
	glVertex3f(  0.0,  0.0,  1.1);
	glVertex3f(  0.1, -0.1,  0.9);
	glVertex3f(  0.0,  0.0,  1.1);
	glVertex3f( -0.1, -0.1,  0.9);

	// left arrow
	glVertex3f(  0.0,  0.0, -1.1);
	glVertex3f(  0.1,  0.1, -0.9);
	glVertex3f(  0.0,  0.0, -1.1);
	glVertex3f( -0.1,  0.1, -0.9);
	glVertex3f(  0.0,  0.0, -1.1);
	glVertex3f(  0.1, -0.1, -0.9);
	glVertex3f(  0.0,  0.0, -1.1);
	glVertex3f( -0.1, -0.1, -0.9);
	glEnd ();

	// right arrow
	glColor4f(std::min(1.0f,r+0.2f), std::min(1.0f,g+0.2f), std::min(1.0f,b+0.2f), 0.5);
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0,  1.1);
	glVertex3f(  0.1,  0.1,  0.9);
	glVertex3f( -0.1,  0.1,  0.9);
	glVertex3f( -0.1, -0.1,  0.9);
	glVertex3f(  0.1, -0.1,  0.9);
	glVertex3f(  0.1,  0.1,  0.9);
	glEnd();
	// left arrow
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f(  0.0,  0.0, -1.1);
	glVertex3f(  0.1,  0.1, -0.9);
	glVertex3f( -0.1,  0.1, -0.9);
	glVertex3f( -0.1, -0.1, -0.9);
	glVertex3f(  0.1, -0.1, -0.9);
	glVertex3f(  0.1,  0.1, -0.9);
	glEnd();
}

void ForceEditPlugin::DrawCircle(float r, float g, float b)
{
	int nside =32;
	const double pi2 = 3.14159265 * 2.0;

	glColor4f(r,g,b,1.0);
	glBegin (GL_LINE_LOOP);
	for (double i = 0; i < nside; i++) 
	{
		glNormal3d (cos (i * pi2 / nside), sin (i * pi2 / nside), 0.0);
		glVertex3d (cos (i * pi2 / nside), sin (i * pi2 / nside), 0.0);
	}
	glEnd();

	glColor4f(std::min(1.0f,r+0.2f), std::min(1.0f,g+0.2f), std::min(1.0f,b+0.2f), 0.5);
	glBegin (GL_TRIANGLE_FAN);
	glVertex3d (0.0, 0.0, 0.0);
	int renderangle;
	if (inputValue>=0)
		renderangle = int(inputValue)%360;
	else
		renderangle = 360 - (int(-inputValue)%360);

	for (double i = 0; i<=renderangle; i++) 
	{
		glVertex3d (cos (i * pi2 / 360.0), sin (i * pi2 / 360.0), 0.0);
	}
	glEnd();
}
void ForceEditPlugin::DrawMeshBox(MeshModel &model)
{
	Box3f b;
	b = model.cm.bbox;
	Point3f mi=b.min;
	Point3f ma=b.max;
	Point3f d3=(b.max-b.min)/4.0;
	Point3f zz(0,0,0);

	// setup
	glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POINT_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(1.0);
	glColor(Color4b::Cyan);

	glPushMatrix();
	glMultMatrix(model.cm.Tr);

	glBegin(GL_LINES);
	glColor3f(1.0, 0.5, 0.5); glVertex3f(mi[0],mi[1],mi[2]); glVertex3f(mi[0]+d3[0],mi[1]+zz[1],mi[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(mi[0],mi[1],mi[2]); glVertex3f(mi[0]+zz[0],mi[1]+d3[1],mi[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(mi[0],mi[1],mi[2]); glVertex3f(mi[0]+zz[0],mi[1]+zz[1],mi[2]+d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(ma[0],mi[1],mi[2]); glVertex3f(ma[0]-d3[0],mi[1]+zz[1],mi[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(ma[0],mi[1],mi[2]); glVertex3f(ma[0]+zz[0],mi[1]+d3[1],mi[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(ma[0],mi[1],mi[2]); glVertex3f(ma[0]+zz[0],mi[1]+zz[1],mi[2]+d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(mi[0],ma[1],mi[2]); glVertex3f(mi[0]+d3[0],ma[1]+zz[1],mi[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(mi[0],ma[1],mi[2]); glVertex3f(mi[0]+zz[0],ma[1]-d3[1],mi[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(mi[0],ma[1],mi[2]); glVertex3f(mi[0]+zz[0],ma[1]+zz[1],mi[2]+d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(ma[0],ma[1],mi[2]); glVertex3f(ma[0]-d3[0],ma[1]+zz[1],mi[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(ma[0],ma[1],mi[2]); glVertex3f(ma[0]+zz[0],ma[1]-d3[1],mi[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(ma[0],ma[1],mi[2]); glVertex3f(ma[0]+zz[0],ma[1]+zz[1],mi[2]+d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(mi[0],mi[1],ma[2]); glVertex3f(mi[0]+d3[0],mi[1]+zz[1],ma[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(mi[0],mi[1],ma[2]); glVertex3f(mi[0]+zz[0],mi[1]+d3[1],ma[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(mi[0],mi[1],ma[2]); glVertex3f(mi[0]+zz[0],mi[1]+zz[1],ma[2]-d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(ma[0],mi[1],ma[2]); glVertex3f(ma[0]-d3[0],mi[1]+zz[1],ma[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(ma[0],mi[1],ma[2]); glVertex3f(ma[0]+zz[0],mi[1]+d3[1],ma[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(ma[0],mi[1],ma[2]); glVertex3f(ma[0]+zz[0],mi[1]+zz[1],ma[2]-d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(mi[0],ma[1],ma[2]); glVertex3f(mi[0]+d3[0],ma[1]+zz[1],ma[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(mi[0],ma[1],ma[2]); glVertex3f(mi[0]+zz[0],ma[1]-d3[1],ma[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(mi[0],ma[1],ma[2]); glVertex3f(mi[0]+zz[0],ma[1]+zz[1],ma[2]-d3[2]);

	glColor3f(1.0, 0.5, 0.5); glVertex3f(ma[0],ma[1],ma[2]); glVertex3f(ma[0]-d3[0],ma[1]+zz[1],ma[2]+zz[2]);
	glColor3f(0.5, 1.0, 0.5); glVertex3f(ma[0],ma[1],ma[2]); glVertex3f(ma[0]+zz[0],ma[1]-d3[1],ma[2]+zz[2]);
	glColor3f(0.5, 0.5, 1.0); glVertex3f(ma[0],ma[1],ma[2]); glVertex3f(ma[0]+zz[0],ma[1]+zz[1],ma[2]-d3[2]);
	glEnd();

	// restore
	glPopMatrix();
	glPopAttrib();
}
//---------------------------------------------------------------------------------------
void ForceEditPlugin::DrawTranslateManipulators(ToothNode* node,GLArea* gla)
{  
	MeshModel& model = *(node->m);
	glPushMatrix();

	Point3f mesh_boxcenter, mesh_origin, mesh_xaxis, mesh_yaxis, mesh_zaxis, new_mesh_origin;
	mesh_boxcenter = original_Transform * model.cm.bbox.Center();
	mesh_origin = original_Transform.GetColumn3(3);
	new_mesh_origin = model.cm.Tr.GetColumn3(3);
	mesh_xaxis = original_Transform.GetColumn3(0);
	mesh_yaxis = original_Transform.GetColumn3(1);
	mesh_zaxis = original_Transform.GetColumn3(2);
	float manipsize = model.cm.bbox.Diag() / 2.0;
	Matrix44f track_rotation;
	gla->trackball.track.rot.ToMatrix(track_rotation);

	glLineWidth(2.0);
	
	glTranslate(model.cm.bbox.Center());
	
	Point3f z(0,0,1);
	Point3f az = z ^ node->pca[0];
	float anglez = acos(z * node->pca[0]) * 180 / FPI;

	Point3f ay = z ^ node->pca[2];
	float angley = acos(z * node->pca[2]) * 180 / FPI;

	Point3f ax = z ^ node->pca[1];
	float anglex = acos(z * node->pca[1]) * 180 / FPI;

	switch(curMode) 
	{
	case ManuMode::None:
		glTranslate(new_mesh_origin);      
		glScale(manipsize);
		glMultMatrix(Inverse(track_rotation));
		glRotatef (90, 0, 1, 0);
		DrawArrows(1.0,0.8,0.5);
		glRotatef (90, 1, 0, 0);
		DrawArrows(1.0,0.8,0.5);
		break;
	case ManuMode::TranX:
		glTranslate(new_mesh_origin);
		glScale(manipsize);
		glRotatef (anglex, ax.X(), ax.Y(), ax.Z());
		DrawArrows(1.0,0,0);
		break;
	case ManuMode::TranY:
		glTranslate(new_mesh_origin);
		glScale(manipsize);
		glRotatef (angley, ay.X(), ay.Y(), ay.Z());
		DrawArrows(0,1.0,0);
		break;
	case ManuMode::TranZ:
		glTranslate(new_mesh_origin);
		glScale(manipsize);
		glRotatef (anglez, az.X(), az.Y(), az.Z());
		DrawArrows(0,0,1.0);
		break;
	}

	glLineWidth(1.0);
	glPopMatrix();
}

void ForceEditPlugin::DrawRotateManipulators(ToothNode* node, GLArea* gla)
{ 
	MeshModel& model = *(node->m);
	glPushMatrix();

	Point3f mesh_boxcenter, mesh_origin, mesh_xaxis, mesh_yaxis, mesh_zaxis, new_mesh_origin;
	mesh_boxcenter = original_Transform * model.cm.bbox.Center();
	mesh_origin = original_Transform.GetColumn3(3);
	new_mesh_origin = model.cm.Tr.GetColumn3(3);
	mesh_xaxis = original_Transform.GetColumn3(0);
	mesh_yaxis = original_Transform.GetColumn3(1);
	mesh_zaxis = original_Transform.GetColumn3(2);
	float manipsize = model.cm.bbox.Diag() / 2.0;
	Matrix44f track_rotation;
	gla->trackball.track.rot.ToMatrix(track_rotation);

	glLineWidth(2.0);
	glTranslate(model.cm.bbox.Center());

	Point3f z(0,0,1);
	Point3f az = z ^ node->pca[0];
	float anglez = acos(z * node->pca[0]) * 180 / FPI;

	Point3f ay = z ^ node->pca[2];
	float angley = acos(z * node->pca[2]) * 180 / FPI;

	Point3f ax = z ^ node->pca[1];
	float anglex = acos(z * node->pca[1]) * 180 / FPI;

	switch(curMode) 
	{
	case ManuMode::None:
		glTranslate(mesh_origin);
		glScale(manipsize);
		glMultMatrix(Inverse(track_rotation));
		DrawCircle(1.0,0.8,0.5);
		break;
	case ManuMode::RotX:
		glTranslate(mesh_origin);
		glScale(manipsize);
		glRotatef (anglex, ax.X(), ax.Y(), ax.Z());
		DrawCircle(1.0,0,0);
		break;
	case ManuMode::RotY:
		glTranslate(mesh_origin);
		glScale(manipsize);
		glRotatef (angley, ay.X(), ay.Y(), ay.Z());
		DrawCircle(0,1.0,0);
		break;
	case ManuMode::RotZ:
		glTranslate(mesh_origin);
		glScale(manipsize);
		glRotatef (anglez, az.X(), az.Y(), az.Z());
		DrawCircle(0,0,1.0);
		break;
	}

	glLineWidth(1.0);
	glPopMatrix();
}

void ForceEditPlugin::DrawManipulators(ToothNode* node,GLArea* gla)
{
	// setup
	glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POINT_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);

	if (curMode == ManuMode::TranX || curMode == ManuMode::TranY || curMode == ManuMode::TranZ)
	{
		DrawTranslateManipulators(node,gla);
	}
	else if (curMode == ManuMode::RotX||curMode==ManuMode::RotY||curMode==ManuMode::RotZ)
	{
		DrawRotateManipulators(node,gla);
	}

	// restore
	glPopAttrib();
}

//---------------------------------------------------------------------------------------
bool ForceEditPlugin::MyPick(const int &x, const int &y, Point3f &pp, float mydepth)
{
	double res[3];
	GLdouble mm[16],pm[16]; GLint vp[4];
	glGetDoublev(GL_MODELVIEW_MATRIX,mm);
	glGetDoublev(GL_PROJECTION_MATRIX,pm);
	glGetIntegerv(GL_VIEWPORT,vp);

	gluUnProject(x,y,mydepth,mm,pm,vp,&res[0],&res[1],&res[2]);
	pp=Point3f(res[0],res[1],res[2]);
	return true;
}

void ForceEditPlugin::Decorate(MeshModel& m,GLArea* gla,QPainter* p)
{ 
	/*
	if (showP.size()>0)
	{
		glPushMatrix();
		
		glPointSize(4.0);
		glBegin(GL_POINTS);
		glColor4f(0,255,0,255);
		glVertex(showP[0]);
		for (int i = 1; i<showP.size()-1;i++)
		{
			glColor4f(255,0,0,255);
			glVertex(showP[i]);
		}
		glColor4f(0,0,255,255);
		glVertex(showP[showP.size()-1]);
		glEnd();
		glPopMatrix();
	}
	*/
	
	if (CSs.size()>0)
	{
		glPushMatrix();
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPointSize(10.0);
		glBegin(GL_POINTS);
		glColor4f(0,255,0,255);
		for (int i = 0; i<CSs.size();i++)
		{
			glColor4f(255,0,0,255);
			glVertex(CSs[i]);
		}
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}

	Point3f center, right, top, front;

	MyPick(gla->width()*0.5, gla->height()*0.5, center, 0.5);
	MyPick(gla->width()*0.99, gla->height()*0.5, right, 0.5);
	MyPick(gla->width()*0.5, gla->height()*0.01, top, 0.5);
	MyPick(gla->width()*0.5, gla->height()*0.5, front, 0.01);

	screen_xaxis = (right - center) * 2.0;
	screen_yaxis = (top - center)   * 2.0;
	screen_zaxis = (front - center) * 2.0;

	if(curMode != ManuMode::None)
	{
		if(curMode != ManuMode::None)
		{
			forceui->ui.manu_value->setValue(displayOffset);
		}
	}

	//选择
	if(CMN)
	{
		DrawMeshBox(*(CMN->m));
		glPushMatrix();
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPointSize(10.0);

		for (int i = 0; i < CMN->frontControlPOrign.size();i++)
		{
			glBegin(GL_POINTS);
			glColor4f(0,255,0,255);
			glVertex(CMN->frontControlPOrign[i]);
			glEnd();
		}
		for (int i = 0; i < CMN->backControlPOrigin.size();i++)
		{
			glBegin(GL_POINTS);
			glColor4f(255,0,0,255);
			glVertex(CMN->backControlPOrigin[i]);
			glEnd();
		}
		glPopAttrib();
		glPopMatrix();
		if (CMN->isattachment)
		{
			AttachmentNode *node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
			forceui->ui.spn_startStep->setValue(node->duration[0]);
			forceui->ui.spn_endStep->setValue(node->duration[1]);
		}
	}

	//新建附件
	if(!CMN->isattachment&&pickloc)
	{
		CFaceO* face1;
		bool result = GLPickTri<CMeshO>::PickNearestFace(curpos.x(),gla->height()-curpos.y(),CMN->m->cm, face1);
		if(result)
		{
			BuildAttach(face1); 
		}
		else
		{
			//QMessageBox::warning(NULL, "warning", "No pickloc,Retry Please!", QMessageBox::Yes, QMessageBox::Yes);
		}
		pickloc = false;
	}

	//附件的平移和旋转
	if(CMN->isattachment)
	{
		DrawManipulators(CMN,gla);
	}

	//右键修改
	if(CMN->isattachment&&isRightClicked)
	{
		popAttAdjMenu(curRightPos);
		isRightClicked = false;
	}
 
	//显示力矩
	if(hasnewfata&&!forceui->ui.cb_hideTorque->isChecked())
	{
		ToothNode* tni;
		if(CMN->isattachment)
		{
			AttachmentNode* ani = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
			tni = ani->parent;
		}
		else
			tni = CMN;
		if(tni->n_TargetForce!=vcg::Point3f(0,0,0))
			DrawVector(tni->n_TargetForce.Normalize()*5,vcg::Color4b::LightGreen);
		if(tni->n_TargetTorque!=vcg::Point3f(0,0,0))
			DrawVector(tni->n_TargetTorque.Normalize()*5,vcg::Color4b::Yellow);
		for(int i = 0;i< tni->attforces.size();i++){
			DrawWireCone(tni->attloc[i],tni->attforces[i].Normalize());
		}
	}

	//测量
	if(is_measure)
	{
		rubberband.Render(this->gla);
		if(rubberband.IsReady())
		{
			Point3f a,b;
			rubberband.GetPoints(a,b);
			vcg::glLabel::render(p,b,QString("%1").arg(Distance(a,b)));
		}
		was_ready=true;
	}

	//if (forceJudge&&CMN->move_path.size()>0)
	//{
	//	forceui->ui.XTtext->setText(QString("%1").arg(CMN->AlignerTorque.X()));
	//	forceui->ui.YTtext->setText(QString("%1").arg(CMN->AlignerTorque.Y()));
	//	forceui->ui.ZTtext->setText(QString("%1").arg(CMN->AlignerTorque.Z()));

	//	forceui->ui.XFtext->setText(QString("%1").arg(CMN->AlignerForce.X()));
	//	forceui->ui.YFtext->setText(QString("%1").arg(CMN->AlignerForce.Y()));
	//	forceui->ui.ZFtext->setText(QString("%1").arg(CMN->AlignerForce.Z()));
	//	//
	//	forceui->ui.TargetTX->setText(QString("%1").arg(CMN->C_TargetTorque.X()));
	//	forceui->ui.TargetTY->setText(QString("%1").arg(CMN->C_TargetTorque.Y()));
	//	forceui->ui.TargetTZ->setText(QString("%1").arg(CMN->C_TargetTorque.Z()));

	//	forceui->ui.TargetFX->setText(QString("%1").arg(CMN->C_TargetForce.X()));
	//	forceui->ui.TargetFY->setText(QString("%1").arg(CMN->C_TargetForce.Y()));
	//	forceui->ui.TargetFZ->setText(QString("%1").arg(CMN->C_TargetForce.Z()));
	//	//
	//	forceui->ui.TXline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->tx()-CMN->move_path[teeth_tree[current_teeth]->current_step]->tx()));
	//	forceui->ui.TXline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->ty()-CMN->move_path[teeth_tree[current_teeth]->current_step]->ty()));
	//	forceui->ui.TXline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->tz()-CMN->move_path[teeth_tree[current_teeth]->current_step]->tz()));
	//	forceui->ui.RXline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->rx()-CMN->move_path[teeth_tree[current_teeth]->current_step]->rx()));
	//	forceui->ui.RYline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->ry()-CMN->move_path[teeth_tree[current_teeth]->current_step]->ry()));
	//	forceui->ui.RZline->setText(QString("%1").arg(CMN->move_path[teeth_tree[current_teeth]->current_step+1]->rz()-CMN->move_path[teeth_tree[current_teeth]->current_step]->rz()));
	//}

	assert(!glGetError());
}

void ForceEditPlugin::UpdateMatrix(MeshModel& model, GLArea* _gla, bool applymouseoffset, bool useinputnumber)
{ 
	Matrix44f newmatrix;

	Matrix44f old_rotation;  
	Matrix44f old_translation;  
	Matrix44f old_meshcenter;
	Matrix44f old_meshuncenter;

	Point3f new_scale;
	Point3f axis;
	float mouseXoff;
	float mouseYoff;

	Point3f mesh_boxcenter, mesh_origin, mesh_xaxis, mesh_yaxis, mesh_zaxis;
	mesh_boxcenter = model.cm.bbox.Center();
	mesh_origin = original_Transform.GetColumn3(3);

	ToothNode* node = teeth_tree[current_teeth]->GetAttByName(model.label());
	mesh_xaxis = node->pca[1];
	mesh_yaxis = node->pca[2];
	mesh_zaxis = node->pca[0];

	AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(model.label());
	if (anode)
	{
		if ((AttType::attRootType(anode->attType) == AttType::optRoot||AttType::attRootType(anode->attType) == AttType::optRot)&&anode->qnodename.length()==5)
		{
			AttachmentNode* pNode = teeth_tree[current_teeth]->GetAttByName(anode->qnodename.left(4));
			if (pNode)
			{
				mesh_boxcenter = pNode->m->cm.bbox.Center();
			}
		}
	}

	delta_Transform.SetIdentity();
	newmatrix.SetIdentity();

	if(curMode == ManuMode::None)
	{
		model.cm.Tr = original_Transform;
	}
	else
	{
		if(curMode != ManuMode::None)  // transform on one axis only
		{
			switch(curMode)          // which axis is active
			{
			case ManuMode::TranX:
				axis = mesh_xaxis;
				break;
			case ManuMode::TranY:
				axis = mesh_yaxis;
				break;
			case ManuMode::TranZ:
				axis = mesh_zaxis;
				break;
			case ManuMode::RotX:
				axis = mesh_xaxis;
				break;
			case ManuMode::RotY:
				axis = mesh_yaxis;
				break;
			case ManuMode::RotZ:
				axis = mesh_zaxis;
				break;
			default: axis = Point3f(1.0, 1.0, 1.0); // it should never arrive here, anyway
			}

			if(curMode == ManuMode::TranX || curMode == ManuMode::TranY || curMode == ManuMode::TranZ)
			{
				// mouse offset -> single axis translation
				float xsign = ((screen_xaxis*axis)>0.0)?1.0:-1.0;
				float ysign = ((screen_yaxis*axis)>0.0)?1.0:-1.0;
				mouseXoff = xsign * screen_xaxis.Norm() * (currScreenOffset_X/float(gla->width()));
				mouseYoff = ysign * screen_yaxis.Norm() * (currScreenOffset_Y/float(gla->height()));
				displayOffset = currOffset + mouseXoff + mouseYoff;

				if(useinputnumber)
					displayOffset = inputValue;

				delta_Transform.SetTranslate(axis * displayOffset);  
				newmatrix = delta_Transform * original_Transform;
			}
			else if(curMode == ManuMode::RotX||curMode==ManuMode::RotY||curMode==ManuMode::RotZ)
			{
				// mouse offset -> single axis rotation
				mouseXoff = (currScreenOffset_X/float(gla->width()));
				mouseYoff = (currScreenOffset_Y/float(gla->height()));
				displayOffset = currOffset + (360.0 * (mouseXoff + mouseYoff));
				if((displayOffset > 360.0) || (displayOffset < -360.0))
					displayOffset = 360.0;
				
				if(useinputnumber)
					displayOffset = inputValue;

				delta_Transform.SetRotateDeg(displayOffset, axis);

				old_rotation = original_Transform;
				old_rotation.SetColumn(3, Point3f(0.0, 0.0, 0.0));
				old_translation.SetTranslate(original_Transform.GetColumn3(3));
				old_meshcenter.SetTranslate(old_rotation * (-mesh_boxcenter));
				old_meshuncenter.SetTranslate(old_rotation * mesh_boxcenter);
				newmatrix = old_translation * old_meshuncenter * delta_Transform * old_meshcenter * old_rotation;
			}
			else
				newmatrix = original_Transform;  // it should never arrive here, anyway
		}
		model.cm.Tr = newmatrix;
	}

	if(applymouseoffset)
	{
		// user finished dragging... accumulation of mouse offset into current offset
		currOffset = displayOffset;
		currOffset_X = displayOffset_X;
		currOffset_Y = displayOffset_Y;
		currOffset_Z = displayOffset_Z;
	}
}
//-----------------------------------------------------------------------------------------
void ForceEditPlugin::mousePressEvent(QMouseEvent* event, MeshModel&, GLArea* gla )
{
	if (event->button() == Qt::RightButton)
	{
		curRightPos = event->globalPos();
		isRightClicked = true;
	}

	if(is_measure)
	{
		if(was_ready||rubberband.IsReady())
		{
			if(firstclick)
			{
				rubberband.Reset();
				firstclick = false;
			}
			else
				firstclick = true;
			was_ready = false;
		}
	}


	//黏贴附件的位置
	if(pickloc)
	{
		curpos = event->pos();
	}
	isMoving = true;
	startdrag = Point2i(event->x(),event->y());

	gla->update();
}

void ForceEditPlugin::mouseMoveEvent(QMouseEvent* event, MeshModel& m, GLArea* gla)
{
	if(is_measure)
		rubberband.Drag(event->pos());

	if(isMoving)
	{
		enddrag = Point2i(event->x(),event->y());
		currScreenOffset_X = enddrag[0] - startdrag[0];
		currScreenOffset_Y = enddrag[1] - startdrag[1];
		if (curMode!=ManuMode::None)
		{
			AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(m.label());
			if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
			{
				QString aname="";
				if (m.label().length()==4)
				{
					aname = m.label()+"0";
					AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
					if (anode)UpdateMatrix(*(anode->m),gla,false);
					UpdateMatrix(m, gla, false);
					
				}
				else if (m.label().length()==5)
				{
					aname = m.label().left(4);
					AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
					if (anode)UpdateMatrix(*(anode->m),gla,false);
					UpdateMatrix(m, gla, false);	
				}
			}
			else
				UpdateMatrix(m, gla, false);
		}	
	}
	gla->update();
}

void ForceEditPlugin::mouseReleaseEvent(QMouseEvent* event, MeshModel& m, GLArea* gla)
{
	if(is_measure)
	{
		rubberband.Pin(event->pos());
	}

	if(isMoving)
	{
		isMoving = false;
		enddrag = Point2i(event->x(),event->y());
		currScreenOffset_X = enddrag[0] - startdrag[0];
		currScreenOffset_Y = enddrag[1] - startdrag[1];
		if (curMode!=ManuMode::None)
		{
			UpdateMatrix(m, gla, false);
		}
	}
	gla->update();
}
//-----------------------------------------------------------------------------------------
void ForceEditPlugin::deleteAtt(AttachmentNode* node)
{
	this->md->delMesh(node->m);
	this->SetCurrentNode(node->parent);
	for (int i = 0; i < teeth_tree[current_teeth]->attNodes.size();i++)
	{
		if (teeth_tree[current_teeth]->attNodes[i]->qnodename == node->qnodename)
		{
			teeth_tree[current_teeth]->attNodes.removeAt(i);
			break;
		}
	}
	this->forceui->rebuildTree();
	this->gla->update();
}

void ForceEditPlugin::slotAdjustParam(int flag)
{
	AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
	if (node == NULL)return;

	vcg::Point3f mp;
	if (node->mp==vcg::Point3f(0,0,0))
	{
		mp = node->barycentric_coord;
	}
	else
		mp= node->mp;

	vcg::Point3f nv;
	if (node->nv == vcg::Point3f(0,0,0))
	{
		nv = node->pca[0];
	}
	else
		nv = node->nv;

	int attType = node->attType;
	ToothNode* parent = node->parent;
	int duration[2];
	duration[0] = node->duration[0];
	duration[1] = node->duration[1];
	QString name = node->qnodename;
	int id = node->id;

	bool isAutoAtt = false;
	Point3f vl = Point3f(0,0,0);
	Point3f vvx = Point3f(0,0,0);
	Point3f vn = Point3f(0,0,0);
	//Point3f p = Point3f(0,0,0);
	Point3f p = mp;
	if(node->isAuto)
	{
		isAutoAtt = node->isAuto;
		vl = node->vl;
		vvx = node->vvx;
		vn = node->vn;
		p = node->p;
	}

	deleteAtt(node);
	MeshModel* m = md->addNewMesh(name,false);
	AttachmentNode *newNode = new AttachmentNode(m,id,name,parent);
	newNode->setType(attType);
	newNode->setPickedInfo(mp,nv);
	newNode->setDuration(duration[0],duration[1]);
	if (isAutoAtt)
	{
		newNode->genAttMesh(showP,flag,0,vl);
	}
	else
		newNode->genAttMesh(showP,flag);
	if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm)
		if (isAutoAtt)
		{
			newNode->initPosition(p,vn,vvx);
		}
		else
			newNode->initPosition(mp,nv);
	if (isAutoAtt)
	{
		newNode->updateColor(vcg::Point3f(0,0,0));
	}
	else
		newNode->updateColor(nv);
	newNode->updateRotation();
	if (isAutoAtt)
	{
		newNode->isAuto = true;
		newNode->p = p;
		newNode->vl = vl;
		newNode->vvx = vvx;
		newNode->vn = vn;
	}
	teeth_tree[current_teeth]->attNodes.push_back(newNode);


	if (AttType::attRootType(newNode->attType) == AttType::optRot)
	{
		AttachmentNode* enode = teeth_tree[current_teeth]->GetAttByName(name+"0");
		if (enode) deleteAtt(enode);

		MeshModel* mm = md->addNewMesh(name+"0",false);
		AttachmentNode* newNNode = new AttachmentNode(mm,teeth_tree[current_teeth]->attNodes.size(),name+"0",parent);
		newNNode->setType(attype);
		newNNode->setPickedInfo(mp,nv);
		newNNode->setDuration(duration[0],duration[1]);
		if (isAutoAtt)
		{
			newNNode->genAttMesh(showP,0,1,vl);
		}
		else
			newNNode->genAttMesh(showP,0,1);
		if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm)
			if (isAutoAtt)
			{
				newNNode->initPosition(p,vn,vvx);
			}
			else
				newNNode->initPosition(mp,nv);
		if (isAutoAtt)
		{
			newNNode->updateColor(vcg::Point3f(0,0,0));
		}
		else
			newNNode->updateColor(nv);
		newNNode->updateRotation();
		if (isAutoAtt)
		{
			newNNode->isAuto = true;
			newNNode->p = p;
			newNNode->vl = vl;
			newNNode->vvx = vvx;
			newNNode->vn = vn;
		}
		teeth_tree[current_teeth]->attNodes.push_back(newNNode);
	}

	SetCurrentNode(newNode);
	forceui->rebuildTree();
	gla->update();
	return;
}

void ForceEditPlugin::popAttAdjMenu(QPoint pos)
{
	QMenu *menu = new QMenu(this->gla);
	QSignalMapper* mapper = new QSignalMapper;
	
	AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
	int aType = node->attType;
	
	int number = AttType::attRootType(aType) * 10;
	if (AttType::attRootType(aType) == AttType::rect)
	{
		QAction* ac_IncL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加长度0.5")));
		QAction* ac_DecL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少长度0.5")));
		connect(ac_IncL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncL,number);
		mapper->setMapping(ac_DecL,-number);
		QAction* ac_IncWL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加长边宽度0.1")));
		QAction* ac_DecWL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少长边宽度0.1")));
		connect(ac_IncWL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecWL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncWL,number+1);
		mapper->setMapping(ac_DecWL,-(number+1));
		QAction* ac_IncWS = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加短边宽度0.1")));
		QAction* ac_DecWS = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少短边宽度0.1")));
		connect(ac_IncWS,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecWS,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncWS,number+2);
		mapper->setMapping(ac_DecWS,-(number+2));
		QAction* ac_IncH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加高度0.1")));
		QAction* ac_DecH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少高度0.1")));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,number+3);
		mapper->setMapping(ac_DecH,-(number+3));
	}
	if (AttType::attRootType(aType) == AttType::ellip)
	{
		QAction* ac_IncH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加高度0.1")));
		QAction* ac_DecH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少高度0.1")));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,number);
		mapper->setMapping(ac_DecH,-number);
	}
	if (AttType::attRootType(aType) == AttType::optRot)
	{
		QAction* ac_IncA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加角度15")));
		QAction* ac_DecA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少角度15")));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
	}
	if (AttType::attRootType(aType) == AttType::optExt)
	{
		QAction* ac_IncL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加长度1")));
		QAction* ac_DecL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少长度1")));
		connect(ac_IncL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncL,number);
		mapper->setMapping(ac_DecL,-number);
	}
	if (AttType::attRootType(aType) == AttType::powerPoint)
	{
		QAction* ac_IncH = menu->addAction(tr("Inc Height by 0.01mm"));
		QAction* ac_DecH = menu->addAction(tr("Dec Height by 0.01mm"));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,number);
		mapper->setMapping(ac_DecH,-number);
		QAction* ac_IncLR = menu->addAction(tr("Inc L-Radius by 0.01mm"));
		QAction* ac_DecLR = menu->addAction(tr("Dec L-Radius by 0.01mm"));
		connect(ac_IncLR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecLR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncLR,number+1);
		mapper->setMapping(ac_DecLR,-(number+1));
		QAction* ac_IncSR = menu->addAction(tr("Inc S-Radius by 0.01mm"));
		QAction* ac_DecSR = menu->addAction(tr("Dec S-Radius by 0.01mm"));
		connect(ac_IncSR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecSR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncSR,number+2);
		mapper->setMapping(ac_DecSR,-(number+2));
		QAction* ac_IncHR = menu->addAction(tr("Inc H-Radius by 0.01mm"));
		QAction* ac_DecHR = menu->addAction(tr("Dec H-Radius by 0.01mm"));
		connect(ac_IncHR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecHR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncHR,number+3);
		mapper->setMapping(ac_DecHR,-(number+3));
	}
	if (AttType::attRootType(aType) == AttType::powerRidge)
	{
		QAction* ac_IncL = menu->addAction(tr("Inc Length by 0.5mm"));
		QAction* ac_DecL = menu->addAction(tr("Dec Length by 0.5mm"));
		connect(ac_IncL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncL,number);
		mapper->setMapping(ac_DecL,-number);
		QAction* ac_IncW = menu->addAction(tr("Inc Width by 0.05mm"));
		QAction* ac_DecW = menu->addAction(tr("Dec Width by 0.05mm"));
		connect(ac_IncW,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecW,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncW,number+1);
		mapper->setMapping(ac_DecW,-(number+1));
	}
	if (AttType::attRootType(aType) == AttType::powerArea)
	{
		QAction* ac_IncA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加A长度0.1")));
		QAction* ac_DecA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少A长度0.1")));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncB = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加B长度0.1")));
		QAction* ac_DecB = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少B长度0.1")));
		connect(ac_IncB,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecB,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncB,number+1);
		mapper->setMapping(ac_DecB,-(number+1));
		QAction* ac_IncH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加高度0.01")));
		QAction* ac_DecH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少高度0.01")));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,number+2);
		mapper->setMapping(ac_DecH,-(number+2));
	}
	if (AttType::attRootType(aType) == AttType::powerArm)
	{
		QAction* ac_IncL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加长度0.5")));
		QAction* ac_DecL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少长度0.5")));
		connect(ac_IncL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncL,number);
		mapper->setMapping(ac_DecL,-number);
		QAction* ac_IncWL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加长边宽度0.1")));
		QAction* ac_DecWL = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少长边宽度0.1")));
		connect(ac_IncWL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecWL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncWL,number+1);
		mapper->setMapping(ac_DecWL,-(number+1));
		QAction* ac_IncWS = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加短边宽度0.1")));
		QAction* ac_DecWS = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少短边宽度0.1")));
		connect(ac_IncWS,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecWS,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncWS,number+2);
		mapper->setMapping(ac_DecWS,-(number+2));
		QAction* ac_IncH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加高度0.1")));
		QAction* ac_DecH = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少高度0.1")));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,number+3);
		mapper->setMapping(ac_DecH,-(number+3));
	}
	if (AttType::attRootType(aType) == AttType::biteRamp)
	{
		QAction* ac_IncL = menu->addAction(tr("Inc Length by 0.2mm"));
		QAction* ac_DecL = menu->addAction(tr("Dec Length by 0.2mm"));
		connect(ac_IncL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecL,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncL,number);
		mapper->setMapping(ac_DecL,-number);
		QAction* ac_IncW = menu->addAction(tr("Inc Width by 0.2mm"));
		QAction* ac_DecW = menu->addAction(tr("Dec Width by 0.2mm"));
		connect(ac_IncW,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecW,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncW,(number+1));
		mapper->setMapping(ac_DecW,-(number+1));
		QAction* ac_IncH = menu->addAction(tr("Inc Height by 0.2mm"));
		QAction* ac_DecH = menu->addAction(tr("Dec Height by 0.2mm"));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,(number+2));
		mapper->setMapping(ac_DecH,-(number+2));
	}
	if(AttType::attRootType(aType) == AttType::buttonCutout)
	{
		QAction* ac_LM = menu->addAction(tr("Left  Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,number);
		mapper->setMapping(ac_RM,-number);
		QAction* ac_IncH = menu->addAction(tr("Inc Height By 0.05mm"));
		QAction* ac_DecH = menu->addAction(tr("Dec Height By 0.05mm"));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,(number+1));
		mapper->setMapping(ac_DecH,-(number+1));
		QAction* ac_IncR = menu->addAction(tr("Inc Radius by 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius by 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+3));
		mapper->setMapping(ac_DecR,-(number+3));
		QAction* ac_IncCR = menu->addAction(tr("Inc circleRadius by 0.01mm"));
		QAction* ac_DecCR = menu->addAction(tr("Dec circleRadius by 0.01mm"));
		connect(ac_IncCR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecCR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncCR,(number+2));
		mapper->setMapping(ac_DecCR,-(number+2));
	}
	if(AttType::attRootType(aType) == AttType::buttonCutout_l)
	{
		QAction* ac_LM = menu->addAction(tr("Left  Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,number);
		mapper->setMapping(ac_RM,-number);
		QAction* ac_IncH = menu->addAction(tr("Inc Height By 0.05mm"));
		QAction* ac_DecH = menu->addAction(tr("Dec Height By 0.05mm"));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,(number+1));
		mapper->setMapping(ac_DecH,-(number+1));
		QAction* ac_IncR = menu->addAction(tr("Inc Radius by 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius by 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+3));
		mapper->setMapping(ac_DecR,-(number+3));
		QAction* ac_IncCR = menu->addAction(tr("Inc circleRadius by 0.01mm"));
		QAction* ac_DecCR = menu->addAction(tr("Dec circleRadius by 0.01mm"));
		connect(ac_IncCR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecCR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncCR,(number+2));
		mapper->setMapping(ac_DecCR,-(number+2));
	}
	if (AttType::attRootType(aType) == AttType::buttonCutout_m)
	{
		QAction* ac_IncH = menu->addAction(tr("Inc horizontal By 0.05mm"));
		QAction* ac_DecH = menu->addAction(tr("Dec horizontal By 0.05mm"));
		connect(ac_IncH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecH,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncH,(number+1));
		mapper->setMapping(ac_DecH,-(number+1));

		QAction* ac_IncV = menu->addAction(tr("Inc vertical By 0.01mm"));
		QAction* ac_DecV = menu->addAction(tr("Dec vertical By 0.01mm"));
		connect(ac_IncV,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecV,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncV,(number+2));
		mapper->setMapping(ac_DecV,-(number+2));
	}
	if (AttType::attRootType(aType) == AttType::leftHook)
	{
		QAction* ac_IncA = menu->addAction(tr("Inc Angle By 5degree"));
		QAction* ac_DecA = menu->addAction(tr("Dec Angle By 5degree"));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncR = menu->addAction(tr("Inc Radius By 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius By 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+1));
		mapper->setMapping(ac_DecR,-(number+1));
		QAction* ac_LM = menu->addAction(tr("Left Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,(number+2));
		mapper->setMapping(ac_RM,-(number+2));

	}
	if (AttType::attRootType(aType) == AttType::rightHook)
	{
		QAction* ac_IncA = menu->addAction(tr("Inc Angle By 5degree"));
		QAction* ac_DecA = menu->addAction(tr("Dec Angle By 5degree"));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncR = menu->addAction(tr("Inc Radius By 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius By 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+1));
		mapper->setMapping(ac_DecR,-(number+1));
		QAction* ac_LM = menu->addAction(tr("Left Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,(number+2));
		mapper->setMapping(ac_RM,-(number+2));
	}
	if (AttType::attRootType(aType) == AttType::leftHook_l)
	{
		QAction* ac_IncA = menu->addAction(tr("Inc Angle By 5degree"));
		QAction* ac_DecA = menu->addAction(tr("Dec Angle By 5degree"));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncR = menu->addAction(tr("Inc Radius By 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius By 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+1));
		mapper->setMapping(ac_DecR,-(number+1));
		QAction* ac_LM = menu->addAction(tr("Left Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,-(number+2));
		mapper->setMapping(ac_RM,(number+2));
	}
	if (AttType::attRootType(aType) == AttType::rightHook_l)
	{
		QAction* ac_IncA = menu->addAction(tr("Inc Angle By 5degree"));
		QAction* ac_DecA = menu->addAction(tr("Dec Angle By 5degree"));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncR = menu->addAction(tr("Inc Radius By 0.01mm"));
		QAction* ac_DecR = menu->addAction(tr("Dec Radius By 0.01mm"));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number+1));
		mapper->setMapping(ac_DecR,-(number+1));
		QAction* ac_LM = menu->addAction(tr("Left Move"));
		QAction* ac_RM = menu->addAction(tr("Right Move"));
		connect(ac_LM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_RM,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_LM,(number+2));
		mapper->setMapping(ac_RM,-(number+2));
	}
	if (AttType::attRootType(aType) == AttType::hook)
	{
		QAction* ac_IncA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("水平增加中心点距离")));
		QAction* ac_DecA = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("水平减少中心点距离")));
		connect(ac_IncA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecA,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncA,number);
		mapper->setMapping(ac_DecA,-number);
		QAction* ac_IncB = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("垂直增加中心点距离")));
		QAction* ac_DecB = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("垂直减少中心点距离")));
		connect(ac_IncB,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecB,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncB,(number+1));
		mapper->setMapping(ac_DecB,-(number+1));
		QAction* ac_IncAn = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加偏转角度")));
		QAction* ac_DecAn = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少偏转角度")));
		connect(ac_IncAn,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecAn,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncAn,(number+2));
		mapper->setMapping(ac_DecAn,-(number+2));
		QAction* ac_IncR = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("增加牵引钩粗度")));
		QAction* ac_DecR = menu->addAction(tr("%1").arg(QString::fromLocal8Bit("减少牵引钩粗度")));
		connect(ac_IncR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		connect(ac_DecR,SIGNAL(triggered(bool)),mapper,SLOT(map()));
		mapper->setMapping(ac_IncR,(number + 3));
		mapper->setMapping(ac_DecR,-(number+ 3));
	}

	connect(mapper,SIGNAL(mapped(int)),this,SLOT(slotAdjustParam(int)));
	menu->exec(pos);
}

void ForceEditPlugin::BuildAttach(CFaceO *face)
 {
	vcg::Point3f mp(0,0,0);
	for(unsigned int j =0;j<3;j++)
	{
		mp += face->V(j)->P();
	}
	mp /= 3;

	int count = 0;
	while (teeth_tree[current_teeth]->GetAttByName(CMN->qnodename+QString::number(count,10)))
	{
		count++;
	}
	QString newNodeName = CMN->qnodename+QString::number(count,10);
	MeshModel* m = md->addNewMesh(newNodeName,false);
	if(attype== AttType::importObj)
	{
		QString  dir, file;
		dir =QCoreApplication::applicationDirPath().toLocal8Bit();

		switch (importID)
		{
		case 0: file = QFileDialog::getOpenFileName(this->gla,tr("Open obj File..."),dir,tr("OBJ (*.obj)"));break;
		case 1: file = dir + "/extraObj/sample/3mm_Rect_Att.obj";											break;
		case 2: file = dir + "/extraObj/sample/4mm_Rect_Att.obj";											break;
		case 3: file = dir + "/extraObj/sample/Opt_Extr_Att_short.obj";										break;
		case 4: file = dir + "/extraObj/sample/Opt_Extr_Att_medium.obj";									break;
		case 5: file = dir + "/extraObj/sample/Opt_Deepbite_Att.obj";										break;
		case 6: file = dir + "/extraObj/sample/Opt_Deepbite&Extrusion_Att.obj";								break;
		case 7: file = dir + "/extraObj/sample/Opt_Rotation_Att_a.obj";										break;
		case 8: file = dir + "/extraObj/sample/Opt_Rotation_Att_b.obj";										break;
		case 9: file = dir + "/extraObj/sample/Opt_Rotation_Att_c.obj";										break;
		case 10:file = dir + "/extraObj/sample/Opt_Rotation_Att_d.obj";										break;
		case 11:file = dir + "/extraObj/sample/Opt_Rotation_Att_e.obj";										break;
		case 12:file = dir + "/extraObj/sample/Opt_Rotation_Att_f.obj";										break;
		case 13:file = dir + "/extraObj/sample/Opt_Root_Control_Att_a.obj";									break;
		case 14:file = dir + "/extraObj/sample/Opt_Root_Control_Att_b.obj";									break;
		case 15:file = dir + "/extraObj/sample/Opt_Root_Control_Att_c1.obj";								break;
		case 16:file = dir + "/extraObj/sample/Opt_Root_Control_Att_c2.obj";								break;
		case 17:file = dir + "/extraObj/sample/G6_Opt_Retraction_Att.obj";									break;
		case 18:file = dir + "/extraObj/sample/G6_Opt_Anchorage_Att.obj";									break;
		case 19:file = dir + "/extraObj/powerarm/3mmCRT.obj";												break;
		case 20:file = dir + "/extraObj/powerarm/4mmCRT.obj";												break;
		case 21:file = dir + "/extraObj/powerarm/5mmCRT.obj";												break;
		case 22:file = dir + "/extraObj/powerarm/3mmCRT_Retention.obj";										break;
		case 23:file = dir + "/extraObj/powerarm/4mmCRT_Retention.obj";										break;
		case 24:file = dir + "/extraObj/powerarm/5mmCRT_Retention.obj";										break;
		case 25:file = dir + "/extraObj/biteramp/Con_BiteRamp_1.obj";										break;
		case 26:file = dir + "/extraObj/biteramp/Con_BiteRamp_2.obj";										break;
		}
		vcg::tri::io::ImporterOBJ<CMeshO>::Info info;
		info.mask = -1;
		vcg::tri::io::ImporterOBJ<CMeshO>::Open (m->cm, file.toAscii().data(), info);
	}
	
	AttachmentNode *newNode = new AttachmentNode(m,teeth_tree[current_teeth]->attNodes.size(),newNodeName,CMN);
	newNode->setType(attype);
	newNode->setPickedInfo(mp,face->N());
	newNode->setDuration(0,teeth_tree[current_teeth]->planstep-1);
	newNode->genAttMesh(showP,0);
	if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm) newNode->initPosition(mp,face->N());
	QString filename="C:\\Users\\User\\Desktop\\73894.obj";
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(newNode->m->cm,filename.toLocal8Bit().data(),0);
	newNode->updateColor(face->N());
	newNode->updateRotation();
	teeth_tree[current_teeth]->attNodes.push_back(newNode);
	
	if (AttType::attRootType(newNode->attType) == AttType::optRot||AttType::attRootType(newNode->attType) == AttType::optRoot)
	{
		MeshModel* mm = md->addNewMesh(newNodeName+"0",false);
		AttachmentNode* newNNode = new AttachmentNode(mm,teeth_tree[current_teeth]->attNodes.size(),newNodeName+"0",CMN);
		newNNode->setType(attype);
		newNNode->setPickedInfo(mp,face->N());
		newNNode->setDuration(0,teeth_tree[current_teeth]->planstep-1);
		newNNode->genAttMesh(showP,0,1);
		if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm)newNNode->initPosition(mp,face->N());
		newNNode->updateColor(face->N());
		newNNode->updateRotation();
		teeth_tree[current_teeth]->attNodes.push_back(newNNode);
	}
	
	SetCurrentNode(newNode);
	forceui->rebuildTree();
	return;
}
void ForceEditPlugin::slotBuildAtt(int _attType)
{
	if (!pickloc)
	{
		pickloc = true;
	}
	attype = _attType;
	if (_attType == AttType::importObj) 
	{
		importID =  forceui->ui.cb_impAtts->currentIndex();
	}
	gla->setCursor(QCursor(QPixmap(":/images/newbox.png"),1,1));
}

//-----------------------------------------------------------------------------------------------
Matrix44f ForceEditPlugin::ArbRot(float t,vcg::Point3f p, vcg::Point3f vv,float alpha)
{
	vcg::Point3f vec=vv.Normalize();
	float u=vec.X();
	float v=vec.Y();
	float w=vec.Z();
	float a=p.X();
	float b=p.Y();
	float c=p.Z();
	Matrix44f ro,ra;
	ro.SetIdentity();
	ra.SetIdentity();
	ro[0][0]=u*u+(v*v+w*w)*cos(t);   ro[0][1]=u*v*(1-cos(t))-w*sin(t); ro[0][2]=u*w*(1-cos(t))+v*sin(t); ro[0][3]=(a*(v*v+w*w)-u*(b*v+c*w))*(1-cos(t))+(b*w-c*v)*sin(t);
	ro[1][0]=u*v*(1-cos(t))+w*sin(t);ro[1][1]=v*v+(u*u+w*w)*cos(t);    ro[1][2]=v*w*(1-cos(t))-u*sin(t); ro[1][3]=(b*(u*u+w*w)-v*(a*u+c*w))*(1-cos(t))+(c*u-a*w)*sin(t);
	ro[2][0]=u*w*(1-cos(t))-v*sin(t);ro[2][1]=v*w*(1-cos(t))+u*sin(t); ro[2][2]=w*w+(u*u+v*v)*cos(t);    ro[2][3]=(c*(u*u+v*v)-w*(a*u+b*v))*(1-cos(t))+(a*v-b*u)*sin(t);

	ra[0][0]=alpha;
	ra[1][1]=alpha;
	ra[2][2]=alpha;
	return ro*ra;
}

vcg::Point3f ForceEditPlugin::NearestVer(ToothNode *tni,vcg::Point3f v)
{
	MeshModel *m=tni->m;
	CMeshO::FaceIterator fi;
	tri::UpdateNormals<CMeshO>::PerVertex(m->cm);
	float dis=100;
	vcg::Point3f N1(0,0,0);
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{  
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();

		if ((v0-v).Norm()+(v1-v).Norm()+(v2-v).Norm()<dis)
		{
			dis=(v0-v).Norm()+(v1-v).Norm()+(v2-v).Norm();
			N1=f->N();
		}

	}
	return N1;
}

vcg::Point3f ForceEditPlugin::NearestVer(ToothNode *tni,Point3f v,Point3f xa){
	//穿过x轴最近的牙齿背面的点
	MeshModel *m=tni->m;
	CMeshO::FaceIterator fi;
	tri::UpdateNormals<CMeshO>::PerVertex(m->cm);
	double dis=100;
	Point3f N1=v;
	/*xa+=v;*/
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)//face.size
	{  
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();
		//double d=sqrt((double)((v0-v).SquaredNorm()-((v0-v)*xa)*((v0-v)*xa)))+sqrt((double)((v1-v).SquaredNorm()-((v1-v)*xa)*((v1-v)*xa)))+sqrt((double)((v2-v).SquaredNorm()-((v2-v)*xa)*((v2-v)*xa)));
		//if (d<dis)
		//{
		//	if ((v0-v)*(v0-v)>0.5)
		// {
		//	 dis=d;
		//	 N1=v0;
		// }
		//}
		vcg::Point3f N=((v0-v1)^(v0-v2)).Normalize();
		float t=(v0-v)*N/(xa*N);
		Point3f P=v+xa*t;
		//判断点在不在三角面片内；
		Point3f E0=v0-P;
		Point3f E1=v1-P;
		Point3f E2=v2-P;
		float s=((E0^E1).Norm()+(E0^E2).Norm()+(E2^E1).Norm())/2.0;
		float s1=((v1-v0)^(v2-v0)).Norm()/2.0;
		if (s==s1 && ((P-v).Norm()>0.1))
		{
			N1=f->N();
			tni->fb=P;
		}
	}
	/*cout<<"点的法向量为"<<N1.X()<<" "<<N1.Y()<<" "<<N1.Z()<<endl;*/
	return N1;

}
void ForceEditPlugin::slotAutoAtt()
{
	foreach(ToothNode* tni,teeth_tree[current_teeth]->nodeList)
	{
		tni->attach.clear();
		for (int j = 0; j < tni->attloc.size(); j++)
		{
			int count = 0;
			while (teeth_tree[current_teeth]->GetAttByName(tni->qnodename + QString::number(count, 10)))
			{
				count++;
			}
			QString newNodeName = tni->qnodename + QString::number(count, 10);
			MeshModel* m = md->addNewMesh(newNodeName, false);
			AttachmentNode *newNode = new AttachmentNode(m, teeth_tree[current_teeth]->attNodes.size(), newNodeName, tni);
			std::string s = tni->nodename;
			Matrix44f ro; ro.SetIdentity();
			newNode->setPickedInfo(tni->attloc[j], tni->attforces[j]);
			Point3f vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1), Ou(0, 0, 0), p = tni->attloc[j], vn(0, 0, 1), vvx(1, 0, 0), vl(0, 0, 0);
			float pra1 = 0, pra2 = 0, alpha = 0.9;
			if (s == "UL6" || s == "UL7" || s == "UR6" || s == "UR7")
			{
				if (tni->toothFeature.crownHeight > 3 && teeth_tree[current_teeth]->hasanchorage)
				{
					if (abs(tni->y_tran) > 2) 
					{
						newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->toothFeature.lateralPivot;
						pra1 = 0; pra2 = 0;
					}
					else if (abs(tni->extint) > 2 || abs(tni->rotation_theta) > 20)
					{
						newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->toothFeature.lateralPivot;
						pra1 = 0; pra2 = 0;
					}
					else if (abs(tni->x_tran) > 2 || abs(tni->x_rotataion) > 20 || abs(tni->y_rotataion) > 20)
					{
						newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
						pra1 = 70; pra2 = 0;
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						vvx = vn^vvx;
					}
					else
					{
						newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot, tni->pca[1]);;
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->fb;
						pra1 = 0; pra2 = 0;
					}
				}
			}
			else if (s == "LL6" || s == "LL7" || s == "LR6" || s == "LR7")
			{
				if (tni->toothFeature.crownHeight > 3 && teeth_tree[current_teeth]->hasanchorage)
				{
					if (abs(tni->y_tran) > 2) 
					{
						newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->toothFeature.lateralPivot;
						pra1 = 0; pra2 = 0;
					}
					else if (abs(tni->extint) > 2 || abs(tni->x_tran) > 2 || abs(tni->y_rotataion) > 20)
					{
						newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->toothFeature.lateralPivot;
						pra1 = 0; pra2 = 0;
					}
					else if (abs(tni->x_rotataion) > 20 || abs(tni->rotation_theta) > 20)
					{
						newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
						pra1 = 70; pra2 = 0;
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						vvx = vn^vvx;
					}
					else {
						newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
						ro = ArbRot(FPi / 2.0, Ou, vy);
						vn = NearestVer(tni, tni->toothFeature.lateralPivot);
						Point3f vv = vx.Normalize() ^ vn.Normalize();
						float rotvalue = acos(vx.Normalize()*vn.Normalize());
						ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
						vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
						p = tni->toothFeature.lateralPivot;
						pra1 = 0; pra2 = 0;
					}
				}
			}
			else if (s == "UL4" || s == "UL5" || s == "UR4" || s == "UR5" || s == "LL4" || s == "LL5" || s == "LR4" || s == "LR5")
			{
				if (abs(tni->extint) > 2 || abs(tni->y_tran) > 2)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->x_tran) > 2)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					vvx = vn^vvx;
				}
				else if (abs(tni->x_rotataion) > 20)
				{
					newNode->setType(AttType::optRoot);//优化控根
					vn = NearestVer(tni, tni->attloc[j]);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv, 0.7);//绕着经过p的vv旋转rotvalue(弧度)

					p = tni->attloc[j];
					pra1 = 0; pra2 = 0;
					if ((p - tni->barycentric_coord)*(tni->pca[2]) < 0)
					{
						vvx = tni->pca[2] - vn / (1.0 / (tni->pca[2] * vn));
					}
					else {
						vvx = -(tni->pca[2] - vn / (1.0 / (tni->pca[2] * vn)));
					}

				}
				else if (abs(tni->rotation_theta) > 20) {
					newNode->setType(AttType::optRotLong);//优化旋转
					vn = NearestVer(tni, tni->attloc[j]);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv, 0.75);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->attloc[j];
					pra1 = 0; pra2 = 0;
					if (tni->attforces[j] * vx < 0)
					{
						if (tni->attforces[j] * vz > 0)
						{
							vl = tni->attforces[j];
						}
						else {
							vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
					else
					{
						if (tni->attforces[j] * vz < 0)
						{
							vl = tni->attforces[j];
						}
						else
						{
							vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
				}
				else if (abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					vvx = vn^vvx;
				}

			}
			else if (s == "UL3" || s == "UR3" || s == "LL3" || s == "LR3")
			{//上下颌尖牙
				if (tni->extint > 2)
				{
					newNode->setType(AttType::optExt);//优化伸长
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (tni->extint < -2)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));//投影到变换后的坐标系下的x坐标
					vvx = vn^vvx;//水平矩形附件再旋转90度
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
				}
				else if (abs(tni->x_tran) > 2 || abs(tni->y_tran) > 2)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->x_rotataion) > 20 || abs(tni->rotation_theta) > 20)
				{
					newNode->setType(AttType::optRotLong);//优化旋转
					vn = NearestVer(tni, tni->attloc[j]);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->attloc[j];
					pra1 = 0; pra2 = 0;
					if (tni->attforces[j] * vx < 0)
					{
						if (tni->attforces[j] * vz > 0)
						{
							vl = tni->attforces[j];
						}
						else {
							vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
					else
					{
						if (tni->attforces[j] * vz < 0)
						{
							vl = tni->attforces[j];
						}
						else
						{
							vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
				}
				else if (abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
					if (tni->attforces[j] * tni->pca[0] < 0)
					{
						vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vn));
					}
					else {
						vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vn)));
					}
				}
			}
			else if (s == "UL1" || s == "UR1" || s == "UL2" || s == "UR2")//上颌切牙
			{
				if (tni->extint > 2)
				{
					newNode->setType(AttType::optExt);//优化伸长
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (tni->extint < -2)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					vvx = vn^vvx;
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
				}
				else if (abs(tni->x_tran) > 2 || abs(tni->y_tran) > 2)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->x_rotataion) > 20)
				{
					newNode->setType(AttType::optRotLong);//优化旋转
					vn = NearestVer(tni, tni->attloc[j]);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->attloc[j];
					pra1 = 0; pra2 = 0;
					if (tni->attforces[j] * vx < 0)
					{
						if (tni->attforces[j] * vz > 0)
						{
							vl = tni->attforces[j];
						}
						else {
							vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
					else
					{
						if (tni->attforces[j] * vz < 0)
						{
							vl = tni->attforces[j];
						}
						else
						{
							vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
							ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
							vl = M(ro, vl);
						}
					}
				}
				else if (abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::powerRidge_b);//powerRidge
					vn=-NearestVer(tni,tni->attloc[j]);
					Point3f vv=vz.Normalize()^vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro=ArbRot(rotvalue,Ou,vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx=tni->pca[0]-vn/(1.0/(tni->pca[0]*vn));
					p=tni->attloc[j];
					pra1=0;pra2=0;
				}
			}
			else if (s == "LL1" || s == "LR1" || s == "LL2" || s == "LR2")//下颌切牙
			{
				if (abs(tni->extint) > 2)
				{
					newNode->setType(AttType::optExt);//优化伸长
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vz.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::powerRidge_b);//powerRidge
					vn=-NearestVer(tni,tni->attloc[j]);
					Point3f vv=vz.Normalize()^vn.Normalize();
					float rotvalue = acos(vz.Normalize()*vn.Normalize());
					ro=ArbRot(rotvalue,Ou,vv);//绕着经过p的vv旋转rotvalue(弧度)
					vvx=tni->pca[0]-vn/(1.0/(tni->pca[0]*vn));
					p=tni->attloc[j];
					pra1=0;pra2=0;
				}
				else
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);//绕y轴旋转90度
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);//找到附着面片的法向量
					Point3f vv = vx.Normalize() ^ vn.Normalize();//计算旋转轴
					float rotvalue = acos(vx.Normalize()*vn.Normalize());//计算选装角度
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
			}
			
			if (newNode->attType<=1) continue;
			newNode->isAuto = true;
			newNode->setDuration(0,teeth_tree[current_teeth]->planstep-1);
			newNode->AttachmentNode::genAttMesh(showP,pra1,0,vl,alpha);
			if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm)
			{ 
				newNode->initPosition(p,vn,vvx);
			}
			newNode->updateColor(vcg::Point3f(0,0,0));
			newNode->updateRotation();
			//newNode->isAuto = true;
			newNode->vl = vl;
			newNode->vvx = vvx;
			newNode->vn = vn;
			newNode->p = p;
			teeth_tree[current_teeth]->attNodes.push_back(newNode);

			if(AttType::attRootType(newNode->attType) == AttType::optRot||AttType::attRootType(newNode->attType) == AttType::optRoot)
			{
				MeshModel *mm = md->addNewMesh(newNodeName+"0",false);
				AttachmentNode* newNNode = new AttachmentNode(mm,teeth_tree[current_teeth]->attNodes.size(),newNodeName + "0",CMN);
				newNNode->setType(newNode->attType);
				newNNode->setPickedInfo(tni->attloc[j], tni->attforces[j]);
				newNNode->setDuration(newNode->duration[0],newNode->duration[1]);
				newNNode->isAuto = true;
				newNNode->genAttMesh(showP,pra1,1,vl,alpha);
				if (AttType::attRootType(newNNode->attType)<AttType::buttonCutout||AttType::attRootType(newNNode->attType)==AttType::biteRamp||AttType::attRootType(newNNode->attType)==AttType::powerPoint)
				{
					newNNode->initPosition(p,vn,vvx);
				}
				newNNode->updateColor(vcg::Point3f(0,0,0));
				newNNode->updateRotation();
				//newNNode->isAuto = true;
				newNNode->vl = vl;
				newNNode->vvx = vvx;
				newNNode->vn = vn;
				newNNode->p = p;
				teeth_tree[current_teeth]->attNodes.push_back(newNNode);
			}
		}
	}
	forceui->rebuildTree();
}
//-------------------------------------------------------------------------------------------------------
void ForceEditPlugin::resetOffsets()
{
	displayOffset = 0;        // mouse offset value (single axis)
	displayOffset_X = 0;      // mouse X offset value
	displayOffset_Y = 0;      // mouse Y offset value
	displayOffset_Z = 0;      // mouse Z offset value
	currOffset = 0;           // combined offset value (single axis)
	currOffset_X = 0;         // X offset value
	currOffset_Y = 0;         // Y offset value
	currOffset_Z = 0;         // Z offset value

	currScreenOffset_X = 0;   // horizontal offset (screen space)
	currScreenOffset_Y = 0;   // vertical offset (screen space)
}

void ForceEditPlugin::FreezeTransfrom()
{
	md->setCurrentMesh(CMN->m->id());
	tri::UpdatePosition<CMeshO>::Matrix(CMN->m->cm, CMN->m->cm.Tr,true);
	tri::UpdateBounding<CMeshO>::Box(CMN->m->cm);
	CMN->m->cm.Tr.SetIdentity();

	AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
	if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
	{
		QString aname="";
		if (CMN->m->label().length()==4)
		{
			aname = CMN->m->label()+"0";
			AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
			if (anode)
			{
				tri::UpdatePosition<CMeshO>::Matrix(anode->m->cm, anode->m->cm.Tr,true);
				tri::UpdateBounding<CMeshO>::Box(anode->m->cm);
				anode->m->cm.Tr.SetIdentity();
			}
		}
		else if (CMN->m->label().length()==5)
		{
			aname = CMN->m->label().left(4);
			AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
			if (anode)
			{
				tri::UpdatePosition<CMeshO>::Matrix(anode->m->cm, anode->m->cm.Tr,true);
				tri::UpdateBounding<CMeshO>::Box(anode->m->cm);
				anode->m->cm.Tr.SetIdentity();
			}
		}
	}

	if (node)
	{
		//printf("[%s %d] %s mp:%f %f %f\n",__FUNCTION__,__LINE__,node->qnodename.toLocal8Bit().data(),
		//	node->mp.X(),node->mp.Y(),node->mp.Z());
		vcg::Point3f baryCenter = vcg::Point3f(0,0,0);
		CMeshO::VertexIterator iter = node->m->cm.vert.begin();
		for (;iter!=node->m->cm.vert.end();iter++)
		{
			baryCenter += (*iter).P();
		}
		baryCenter /= node->m->cm.vert.size();
		node->mp = baryCenter;
		node->p = node->mp;
		//printf("[%s %d] %s mp:%f %f %f\n",__FUNCTION__,__LINE__,node->qnodename.toLocal8Bit().data(),
		//	node->mp.X(),node->mp.Y(),node->mp.Z());
	}


	inputValue=0;
	resetOffsets();

	original_Transform = vcg::Matrix44f::Identity();
	delta_Transform = vcg::Matrix44f::Identity();

	gla->update();
}

void ForceEditPlugin::applyMotion(MeshModel& model, GLArea* _gla)
{
	inputValue=0;
	resetOffsets();
	
	// storing start matrix
	original_Transform = model.cm.Tr;
	delta_Transform = vcg::Matrix44f::Identity();

	gla->update();
}

void ForceEditPlugin::slotManStart()
{
	if (forceui->ui.manu_tranX->isChecked()||forceui->ui.manu_tranY->isChecked()||forceui->ui.manu_tranZ->isChecked()) // translate
	{
		if (forceui->ui.manu_tranX->isChecked()) curMode = ManuMode::TranX;
		if (forceui->ui.manu_tranY->isChecked()) curMode = ManuMode::TranY;
		if (forceui->ui.manu_tranZ->isChecked()) curMode = ManuMode::TranZ;
		resetOffsets();
		AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
		if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
		{
			QString aname="";
			if (CMN->m->label().length()==4)
			{
				aname = CMN->m->label()+"0";
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
			else if (CMN->m->label().length()==5)
			{
				aname = CMN->m->label().left(4);
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
		}
		else
			UpdateMatrix(*(CMN->m),gla,false);
	}
	else if (forceui->ui.manu_rotX->isChecked()||forceui->ui.manu_rotY->isChecked()||forceui->ui.manu_rotZ->isChecked()) // rotate
	{
		if (forceui->ui.manu_rotX->isChecked()) curMode = ManuMode::RotX;
		if (forceui->ui.manu_rotY->isChecked()) curMode = ManuMode::RotY;
		if (forceui->ui.manu_rotZ->isChecked()) curMode = ManuMode::RotZ;
		resetOffsets();
		AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
		if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
		{
			QString aname="";
			if (CMN->m->label().length()==4)
			{
				aname = CMN->m->label()+"0";
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
			else if (CMN->m->label().length()==5)
			{
				aname = CMN->m->label().left(4);
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
		}
		else
			UpdateMatrix(*(CMN->m),gla,false);
	}
	else
	{
		curMode = ManuMode::None;
		resetOffsets();
		AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
		if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
		{
			QString aname="";
			if (CMN->m->label().length()==4)
			{
				aname = CMN->m->label()+"0";
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
			else if (CMN->m->label().length()==5)
			{
				aname = CMN->m->label().left(4);
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)UpdateMatrix(*(anode->m),gla,false);
				UpdateMatrix(*(CMN->m), gla,false);
			}
		}
		else
			UpdateMatrix(*(CMN->m),gla,false);
	}
	gla->update();
}

void ForceEditPlugin::slotManApply()
{
	if(curMode != ManuMode::None)
	{
		AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
		if (node&&(AttType::attRootType(node->attType)==AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
		{
			QString aname="";
			if (CMN->m->label().length()==4)
			{
				aname = CMN->m->label()+"0";
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				if (anode)applyMotion(*(anode->m),gla);
				applyMotion(*(CMN->m), gla);
			}
			else if (CMN->m->label().length()==5)
			{
				aname = CMN->m->label().left(4);
				AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
				applyMotion(*(CMN->m), gla);
				if (anode)applyMotion(*(anode->m),gla);
			}
		}
		else
			applyMotion(*(CMN->m), gla);
	}
	
	FreezeTransfrom();


	curMode = ManuMode::None;
	forceui->ui.manu_none->setChecked(true);
	gla->update();
}

void ForceEditPlugin::slotGetManuValue(double delta)
{
	inputValue = delta;
	AttachmentNode* node = teeth_tree[current_teeth]->GetAttByName(CMN->qnodename);
	if (node&&(AttType::attRootType(node->attType) == AttType::optRot||AttType::attRootType(node->attType)==AttType::optRoot))
	{
		QString aname="";
		if (CMN->m->label().length()==4)
		{
			aname = CMN->m->label()+"0";
			AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
			if (anode)UpdateMatrix(*(anode->m),gla,false,true);
			UpdateMatrix(*(CMN->m), gla,false,true);
		}
		else if (CMN->m->label().length()==5)
		{
			aname = CMN->m->label().left(4);
			AttachmentNode* anode = teeth_tree[current_teeth]->GetAttByName(aname);
			if (anode)UpdateMatrix(*(anode->m),gla,false,true);
			UpdateMatrix(*(CMN->m), gla,false,true);
		}
	}
	else
		UpdateMatrix(*(CMN->m), gla, false,true);

	gla->update();
}

void ForceEditPlugin::slotChangeDuration()
{
	if (CMN->isattachment)
	{
		AttachmentNode *node = static_cast<AttachmentNode*>(CMN);
		node->duration[0] = forceui->ui.spn_startStep->value();
		node->duration[1] = forceui->ui.spn_endStep->value();
		qDebug()<<__FUNCTION__<<__LINE__<<" "<<node->duration[0]<<" "<<node->duration[1];
	}
}

//--------------------------------------------------------------------------------------------------------
//按第一行展开计算|A|
double getA(double arcs[L][L],int n)
{
	if(n==1)
	{
		return arcs[0][0];
	}
	double ans = 0;
	double temp[L][L]={0.0};
	int i,j,k;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n-1;j++)
		{
			for(k=0;k<n-1;k++)
			{
				temp[j][k] = arcs[j+1][(k>=i)?k+1:k];
			}
		}
		double t = getA(temp,n-1);
		if(i%2==0)
		{
			ans += arcs[0][i]*t;
		}
		else
		{
			ans -=  arcs[0][i]*t;
		}
	}
	return ans;
}
//计算每一行每一列的每个元素所对应的余子式，组成A*
void  getAStart(double arcs[L][L],int n,double ans[L][L])
{
	if(n==1)
	{
		ans[0][0] = 1;
		return;
	}
	int i,j,k,t;
	double temp[L][L];
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n-1;k++)
			{
				for(t=0;t<n-1;t++)
				{
					temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
				}
			}
			ans[j][i]  =  getA(temp,n-1);
			if((i+j)%2 == 1)
			{
				ans[j][i] = - ans[j][i];
			}
		}
	}
}

//得到给定矩阵src的逆矩阵保存到des中。
bool GetMatrixInverse(int n,double src[L][L],double des[L][L])
{
	double flag=getA(src,n);
	double t[L][L];
	if(flag==0)
	{
		return false;
	}
	else
	{
		getAStart(src,n,t);
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				des[i][j]=t[i][j]/flag;
			}
		}
	}
	return true;
}

inline vcg::Matrix44f CalRotation(ToothNode *mn, int loc,vcg::Point3f center,double theta = 0.0)
{//计算绕a旋转theta角度的旋转矩阵 
	vcg::Matrix44f res,M,tran,ro ,test;
	res.SetIdentity();
	//M = ro*tran,res = M'*Rz(theta)*M
	tran.SetIdentity();
	ro.SetIdentity();
	M.SetIdentity();

	tran[0][3] = -center.X();
	tran[1][3] = -center.Y();
	tran[2][3] = -center.Z();

	ro[0][0] = mn->old_pca[1].X();ro[0][1] = mn->old_pca[1].Y();ro[0][2] = mn->old_pca[1].Z();
	ro[1][0] = mn->old_pca[2].X();ro[1][1] = mn->old_pca[2].Y();ro[1][2] = mn->old_pca[2].Z();
	ro[2][0] = mn->old_pca[0].X();ro[2][1] = mn->old_pca[0].Y();ro[2][2] = mn->old_pca[0].Z();

	M = ro * tran;
	vcg::Matrix44f rz;
	rz.SetIdentity();
	float cosr = cos(theta*FPi/180);
	float sinr = sin(theta*FPi/180);
	if(loc == 0){
		rz[0][0] = cosr;  rz[0][1] = sinr;  rz[0][2] = 0;
		rz[1][0] = -sinr; rz[1][1] = cosr;  rz[1][2] = 0;
		rz[2][0] = 0;     rz[2][1] = 0;     rz[2][2] = 1;
	}
	else if(loc == 2){
		rz[0][0] = cosr;    rz[0][1] = 0;   rz[0][2] = sinr;
		rz[1][0] = 0;       rz[1][1] = 1;   rz[1][2] = 0;
		rz[2][0] = -sinr;   rz[2][1] = 0;   rz[2][2] = cosr;
	}
	else if(loc == 1){
		rz[0][0] = 1;    rz[0][1] = 0;      rz[0][2] = 0;
		rz[1][0] = 0;    rz[1][1] = cosr;   rz[1][2] = sinr;
		rz[2][0] = 0;    rz[2][1] = -sinr;  rz[2][2] = cosr;
	}

	tran[0][3] = -tran[0][3];
	tran[1][3] = -tran[1][3];
	tran[2][3] = -tran[2][3];
	//res = M * rz * M.transpose();
	res = tran * ro.transpose() * rz * M;
	return res;
}

inline vcg::Matrix44f CalRotation(ToothNode *mn,vcg::Point3f center,float theta1 = 0.0,float theta2 = 0.0,float theta3 = 0.0 )
{//计算绕旋转中心center旋转theta角度的旋转矩阵,这个函数绕x,y,z都转动 
	//修改
	vcg::Matrix44f res,M,tran,ro ,test;
	res.SetIdentity();
	//M = ro*tran,res = M'*Rz(theta)*M
	tran.SetIdentity();
	ro.SetIdentity();
	M.SetIdentity();

	tran[0][3] = -center.X();
	tran[1][3] = -center.Y();
	tran[2][3] = -center.Z();

	ro[0][0] = mn->old_pca[1].X();ro[0][1] = mn->old_pca[1].Y();ro[0][2] = mn->old_pca[1].Z();
	ro[1][0] = mn->old_pca[2].X();ro[1][1] = mn->old_pca[2].Y();ro[1][2] = mn->old_pca[2].Z();
	ro[2][0] = mn->old_pca[0].X();ro[2][1] = mn->old_pca[0].Y();ro[2][2] = mn->old_pca[0].Z();

	M = ro * tran;//将点全部移动到旋转中心进行转动
	vcg::Matrix44f rz,rz1,rz2,rz3;//rz是总体的变换矩阵，不是绕z轴的变换矩阵，具体绕哪一个轴，是根据loc的变化来定
	rz.SetIdentity();
	//绕X轴旋转
	float cosr = cos(theta1 * FPi/180);
	float sinr = sin(theta1 * FPi/180);
	rz[0][0] = 1;    rz[0][1] = 0;      rz[0][2] = 0;
	rz[1][0] = 0;    rz[1][1] = cosr;   rz[1][2] = sinr;
	rz[2][0] = 0;    rz[2][1] = -sinr;  rz[2][2] = cosr;
	rz1=rz;
	//绕Y轴旋转
	cosr = cos(theta2*FPi/180);
	sinr = sin(theta2*FPi/180);
	rz[0][0] = cosr;    rz[0][1] = 0;   rz[0][2] = sinr;
	rz[1][0] = 0;       rz[1][1] = 1;   rz[1][2] = 0;
	rz[2][0] = -sinr;   rz[2][1] = 0;   rz[2][2] = cosr;
	rz2=rz;
	//绕Z轴旋转
	cosr = cos(theta3*FPi/180);
	sinr = sin(theta3*FPi/180);
	rz[0][0] = cosr;  rz[0][1] = sinr;  rz[0][2] = 0;
	rz[1][0] = -sinr; rz[1][1] = cosr;  rz[1][2] = 0;
	rz[2][0] = 0;     rz[2][1] = 0;     rz[2][2] = 1;
	rz3=rz;

	//总的旋转矩阵
	rz=rz1*rz2*rz3;

	tran[0][3] = -tran[0][3];
	tran[1][3] = -tran[1][3];
	tran[2][3] = -tran[2][3];

	res = tran * ro.transpose() * rz * M;
	//主成分矩阵是正交矩阵，其逆就是其转置
	return res;
}

inline vcg::Matrix44f CalTran(ToothNode *mn, int loc = 0,float theta = 0.0)
{
	Matrix44f res;
	res.SetIdentity();
	res[0][3] = theta*mn->old_pca[loc].X();
	res[1][3] = theta*mn->old_pca[loc].Y();
	res[2][3] = theta*mn->old_pca[loc].Z();
	return res;
}

inline vcg::Matrix44f CalTran(ToothNode *mn,float theta1 = 0.0,float theta2 = 0.0,float theta3 = 0.0)
{
	Matrix44f res;
	res.SetIdentity();
	//沿着局部坐标X轴平移得到的全局坐标平移变换矩阵，1：X,2:Y,0:Z
	res[0][3] = theta1*mn->old_pca[1].X();
	res[1][3] = theta1*mn->old_pca[1].Y();
	res[2][3] = theta1*mn->old_pca[1].Z();
	//沿着局部坐标Y轴平移得到的全局坐标平移变换矩阵
	res[0][3] += theta2*mn->old_pca[2].X();
	res[1][3] += theta2*mn->old_pca[2].Y();
	res[2][3] += theta2*mn->old_pca[2].Z();
	//沿着局部坐标Z轴平移得到的全局坐标平移变换矩阵
	res[0][3] += theta3*mn->old_pca[0].X();
	res[1][3] += theta3*mn->old_pca[0].Y();
	res[2][3] += theta3*mn->old_pca[0].Z();

	return res;
}

inline double Area(MeshModel *m)
{
	double S=0;
	CMeshO::FaceIterator fi;
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f vloc=vi->P();
		vcg::Point3f vloc1=vi1->P();
		vcg::Point3f vloc2=vi2->P();
		vcg::Point3f V1=vloc1 - vloc;
		vcg::Point3f V2=vloc2 - vloc;
		vcg::Point3f vx = V1^V2;
		double a = vx.Norm();
		S=S+(a/2.0f);
	}
	return S;
}
inline vcg::Point3f Force(MeshModel *m,Matrix44f res){//计算模型移动一定距离后所需要的力
	vcg::Point3f F;
	double si = Area(m);
	//std::cout << si << std::endl;
	F.X()=Kr*si*res[0][3];
	F.Y()=Kr*si*res[1][3];
	F.Z()=Kr*si*res[2][3];
	return F;
}
inline vcg::Point3f Force(ToothNode *tni,Matrix44f res)
{
	vcg::Point3f F;
	double si = Area(tni->m);
	F.X()=Kr*si*res[0][3];
	F.Y()=Kr*si*res[1][3];
	F.Z()=Kr*si*res[2][3];
	return F;
}

inline Point3f Cf(CFaceO *fi)
{
	CFaceO *f = &(*fi);
	CVertexO *vi=f->V(0);
	CVertexO *vi1=f->V(1);
	CVertexO *vi2=f->V(2);
	vcg::Point3f vloc=vi->P();
	vcg::Point3f vloc1=vi1->P();
	vcg::Point3f vloc2=vi2->P();
	return (vloc+vloc1+vloc2)/(3.0f);
}

inline Point3f Cs(MeshModel *m)
{
	float S=Area(m);
	vcg::Point3f K;
	K.X()=0;
	K.Y()=0;
	K.Z()=0;
	CMeshO::FaceIterator fi;
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)//face.size
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f vloc=vi->P();
		vcg::Point3f vloc1=vi1->P();
		vcg::Point3f vloc2=vi2->P();
		vcg::Point3f V1=vloc-vloc1;
		vcg::Point3f V2=vloc-vloc2;
		K.X()+=(V1^V2).Norm()*Cf(f).X()/(2.0f);
		K.Y()+=(V1^V2).Norm()*Cf(f).Y()/(2.0f);
		K.Z()+=(V1^V2).Norm()*Cf(f).Z()/(2.0f);
	}
	return K/S;
}

inline Point3f Cs(ToothNode *TN)//计算模型的中心
{
	vcg::Point3f K;
	K=Cs(TN->m);//获取重心位置
	K-=(TN->pca[0])*(TN->toothFeature.crownHeight)*2.0/3.0;//向局部坐标系的z轴反方向移动2/3的冠高
	CSs.push_back(K);
	return K;

}

inline vcg::Point3f Torque(ToothNode *tni,Matrix44f res)
{//计算旋转了res变换后的力矩
	Point3f T(0,0,0);
	T.SetZero();//初始化
	CMeshO::FaceIterator fi;
	CMeshO::VertexIterator Vi;
	Point3f Csi=Cs(tni);
	for (fi = tni->m->cm.face.begin();fi!=tni->m->cm.face.end();fi++)//face.size
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P()-Csi;//平移到重心位置
		vcg::Point3f v1=vi1->P()-Csi;
		vcg::Point3f v2=vi2->P()-Csi;
		vcg::Point3f V1=v0-v1;
		vcg::Point3f V2=v0-v2;
		float Sf=(V1^V2).Norm()/(2.0f);
		vcg::Point3f X_1=((v0+v1+v2)^(M(res,v0)+M(res,v1)+M(res,v2)))+(v0^M(res,v0))+(v1^M(res,v1))+(v2^M(res,v2));
		T.X()+=Kr/(12.0f)*Sf*X_1.X();
		T.Y()+=Kr/(12.0f)*Sf*X_1.Y();
		T.Z()+=Kr/(12.0f)*Sf*X_1.Z();
	}
	return -T;
}

inline void ChooseAttachFace(MeshModel *m,Point3f *Nm,float Xbound,float Yboundmin,float Yboundmax,float Zboundmin, float Zboundmax,double base[3][3],Point3f TR,bool pickarea) 
{
	//计算附着的面片的区域
	vcg::Color4f c1;
	CMeshO::FaceIterator fi;
	m->updateDataMask(MeshModel::MM_VERTFACETOPO);//获取权限
	m->updateDataMask(MeshModel::MM_FACEFACETOPO);
	m->updateDataMask(MeshModel::MM_FACECOLOR);
	m->updateDataMask(MeshModel::MM_FACEFLAGBORDER);

	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{  
		(*fi).ClearS();//不满足条件就标记
		(*fi).V(0)->ClearS();
		(*fi).V(1)->ClearS();
		(*fi).V(2)->ClearS();
	}

	int flag;
	double Xboundmin=0;
	double Zmax=0;
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();
		v0=M(base,v0-TR);//进行坐标转化
		v1=M(base,v1-TR);
		v2=M(base,v2-TR);
		bool k1=(v0.X()>Xbound)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0);//在边界包围盒内的点
		bool k2=(v1.X()>Xbound)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0);
		bool k3=(v2.X()>Xbound)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0);
		if (k1&&k2&&k3)
		{   
			double Zm=max(max(v0.Z(),v1.Z()),v2.Z());
			double Xb=v0.X()*(v0.Z()==Zm)+v1.X()*(v1.Z()==Zm)+v2.X()*(v2.Z()==Zm);//得到选择好的区域的最高点的X值，使得它成为最精确的边界
			Xboundmin=(Zm>Zmax)?Xb:Xboundmin;
			Zmax=(Zm>Zmax)?Zm:Zmax;

		}
	}
	
	int k=1;
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();
		v0=M(base,v0-TR);//进行坐标转化
		v1=M(base,v1-TR);
		v2=M(base,v2-TR);
		bool k1,k2,k3;
		if (pickarea)
		{
			k1=(v0.X()>Xboundmin)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0)&&((v0.Z()-Zboundmin)*(v0.Z()-Zboundmax)<0);//在边界包围盒内的点
			k2=(v1.X()>Xboundmin)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0)&&((v1.Z()-Zboundmin)*(v1.Z()-Zboundmax)<0);
			k3=(v2.X()>Xboundmin)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0)&&((v2.Z()-Zboundmin)*(v2.Z()-Zboundmax)<0);
		}
		else
		{
			k1=(v0.X()<Xboundmin)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0)&&((v0.Z()-Zboundmin)*(v0.Z()-Zboundmax)<0);//在边界包围盒内的点
			k2=(v1.X()<Xboundmin)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0)&&((v1.Z()-Zboundmin)*(v1.Z()-Zboundmax)<0);
			k3=(v2.X()<Xboundmin)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0)&&((v2.Z()-Zboundmin)*(v2.Z()-Zboundmax)<0);
		}

		if ((k1&&k2&& k3)!=1)
		{    
			(*fi).SetS();//不满足条件就标记
			//(*fi).V(0)->C()= Color4b::Red;
			//(*fi).V(1)->C()= Color4b::Red;
			//(*fi).V(2)->C()= Color4b::Red;
		}else{
			(*fi).C()= Color4b::Red;
			(*fi).V(0)->C()= Color4b::Red;
			(*fi).V(1)->C()= Color4b::Red;
			(*fi).V(2)->C()= Color4b::Red;
			(*fi).V(0)->SetS();
			(*fi).V(1)->SetS();
			(*fi).V(2)->SetS();
			*Nm=*Nm+(*fi).N();
			k++;
		}
	}
	*Nm=*Nm/k;
}
inline vcg::Point3f OptimFT(Point3f d,Point3f T,Point3f F)
{//计算给定的力和力矩最佳的力
	MatrixXd B(3,3);
	B<<
		1+d.Y()*d.Y()+d.Z()*d.Z()+lamda, -d.X()*d.Y() ,-d.X()*d.Z(),
		-d.X()*d.Y(), 1+d.X()*d.X()+d.Z()*d.Z()+lamda, -d.Y()*d.Z(),
		-d.X()*d.Z(), -d.Y()*d.Z(), 1+d.X()*d.X()+d.Y()*d.Y()+lamda;

	MatrixXd b(3,1);
	b<<
		F.X()+d.Z()*T.Y()-d.Y()*T.Z(),
		F.Y()+d.X()*T.Z()-d.Z()*T.X(),
		F.Z()+d.Y()*T.X()-d.X()*T.Y();

	VectorXd Fs(3);
	Fs=B.inverse()*b;
	Point3f Fac(Fs(0),Fs(1),Fs(2));
	return Fac;
}

inline vcg::Point3f OptimFT(Point3f d,Point3f T,Point3f F,Point3f Nm)
{
	MatrixXd B(4,4);
	B<<
		1+d.Y()*d.Y()+d.Z()*d.Z()+lamda, -d.X()*d.Y() ,-d.X()*d.Z(),Nm.X(),
		-d.X()*d.Y(), 1+d.X()*d.X()+d.Z()*d.Z()+lamda, -d.Y()*d.Z(),Nm.Y(),
		-d.X()*d.Z(), -d.Y()*d.Z(), 1+d.X()*d.X()+d.Y()*d.Y()+lamda,Nm.Z(),
		Nm.X(),Nm.Y(),Nm.Z(),0;

	MatrixXd b(4,1);
	b<<
		F.X()+d.Z()*T.Y()-d.Y()*T.Z(),
		F.Y()+d.X()*T.Z()-d.Z()*T.X(),
		F.Z()+d.Y()*T.X()-d.X()*T.Y(),
		0;

	VectorXd Fs(4);
	Fs=B.inverse()*b;
	Point3f Fac(Fs(0),Fs(1),Fs(2));
	return Fac;
}

inline vcg::Point3f OptimFT(Point3f d,Point3f T,Point3f F,Point3f Nm,bool pickarea)
{
	MatrixXd B1(3,3),B2(3,3),B3(3,3),B4(3,3);
	B1<< 1,0,0,
		d.X(),d.Y(),d.Z(),
		Nm.X(),Nm.Y(),Nm.Z();
	B2<< 0,1,0,
		d.X(),d.Y(),d.Z(),
		Nm.X(),Nm.Y(),Nm.Z();
	B3<< 0,0,1,
		d.X(),d.Y(),d.Z(),
		Nm.X(),Nm.Y(),Nm.Z();
	B4<<T.X(),T.Y(),T.Z(),
		d.X(),d.Y(),d.Z(),
		Nm.X(),Nm.Y(),Nm.Z();
	double a=2*(pow(Nm.Norm(),2)+pow(B1.determinant(),2)+pow(B2.determinant(),2)+pow(B3.determinant(),2));
	double b=F*Nm*2+2*B4.determinant()-lamda*pow(Nm.Norm(),2);
	double k=abs(b/a);
	Point3f Fac(0,0,0);
	Fac.X()=k*Nm.X();
	Fac.Y()=k*Nm.Y();
	Fac.Z()=k*Nm.Z();

	return Fac;
}

inline void ChooseAttachFace(ToothNode *node,Point3f *Nm,float Xbound,float Yboundmin,float Yboundmax,float Zboundmin, float Zboundmax,double base[3][3],Point3f TR,bool pickarea)//20170906ÐÞ¸Ä
{
	MeshModel *m=node->m;
	vcg::Color4f c1;
	CMeshO::FaceIterator fi,fi2;
	m->updateDataMask(MeshModel::MM_VERTFACETOPO);
	m->updateDataMask(MeshModel::MM_FACEFACETOPO);
	m->updateDataMask(MeshModel::MM_FACEFLAGBORDER);


	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{  
		(*fi).ClearS();
		(*fi).V(0)->ClearS();
		(*fi).V(1)->ClearS();
		(*fi).V(2)->ClearS();
	}

	int flag;
	double Xboundmin=0;
	double Zmax=0;
	Point3f mc=node->barycentric_coord;
	Point3f x=node->pca[1];
	Point3f y=node->pca[2];
	Point3f z=node->pca[0];
	for (fi = m->cm.face.begin();fi!=m->cm.face.end();fi++)
	{
		CFaceO *f = &(*fi);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();

		v0=M(v0,mc,x,y,z);
		v1=M(v1,mc,x,y,z);
		v2=M(v2,mc,x,y,z);
		bool k1=(v0.X()>Xbound)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0);
		bool k2=(v1.X()>Xbound)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0);
		bool k3=(v2.X()>Xbound)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0);
		if (k1&&k2&&k3)
		{   
			double Zm=max(max(v0.Z(),v1.Z()),v2.Z());
			double Xb=v0.X()*(v0.Z()==Zm)+v1.X()*(v1.Z()==Zm)+v2.X()*(v2.Z()==Zm);
			Xboundmin=(Zm>Zmax)?Xb:Xboundmin;
			Zmax=(Zm>Zmax)?Zm:Zmax;

		}
	}

	int k=1;
	for (fi2 = m->cm.face.begin();fi2!=m->cm.face.end();fi2++)
	{
		CFaceO *f = &(*fi2);
		CVertexO *vi=f->V(0);
		CVertexO *vi1=f->V(1);
		CVertexO *vi2=f->V(2);
		vcg::Point3f v0=vi->P();
		vcg::Point3f v1=vi1->P();
		vcg::Point3f v2=vi2->P();

		v0=M(v0,mc,x,y,z)-TR;
		v0.X()=v0.X()+TR.X();
		v1=M(v1,mc,x,y,z)-TR;
		v1.X()=v1.X()+TR.X();
		v2=M(v2,mc,x,y,z);
		v2.X()=v2.X()+TR.X();
		bool k1,k2,k3;
		if (pickarea)
		{
			k1=(v0.X()>Xboundmin)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0)&&((v0.Z()-Zboundmin)*(v0.Z()-Zboundmax)<0);//ÔÚ±ß½ç°üÎ§ºÐÄÚµÄµã
			k2=(v1.X()>Xboundmin)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0)&&((v1.Z()-Zboundmin)*(v1.Z()-Zboundmax)<0);
			k3=(v2.X()>Xboundmin)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0)&&((v2.Z()-Zboundmin)*(v2.Z()-Zboundmax)<0);
		}
		else
		{   std::string s=node->nodename;
		if (s=="LL1"||s=="LR1"||s=="LL2"||s=="LR2"||s=="UL1"||s=="UR1"||s=="UL2"||s=="UR2")
		{
			Xboundmin=M(node->toothFeature.initlateralPivot,mc,x,y,z).X();
		}
		k1=(v0.X()<Xboundmin)&&((v0.Y()-Yboundmin)*(v0.Y()-Yboundmax)<0)&&((v0.Z()-Zboundmin)*(v0.Z()-Zboundmax)<0);//ÔÚ±ß½ç°üÎ§ºÐÄÚµÄµã
		k2=(v1.X()<Xboundmin)&&((v1.Y()-Yboundmin)*(v1.Y()-Yboundmax)<0)&&((v1.Z()-Zboundmin)*(v1.Z()-Zboundmax)<0);
		k3=(v2.X()<Xboundmin)&&((v2.Y()-Yboundmin)*(v2.Y()-Yboundmax)<0)&&((v2.Z()-Zboundmin)*(v2.Z()-Zboundmax)<0);
		}

		if ((k1&&k2&& k3)!=1)
		{    
			//(*fi2).SetS();
			(*fi2).SetV();
		}
		else
		{
			//(*fi2).V(0)->C()= Color4b::Red;
			//(*fi2).V(1)->C()= Color4b::Red;
			//(*fi2).V(2)->C()= Color4b::Red;
			(*fi2).V(0)->SetS();
			(*fi2).V(1)->SetS();
			(*fi2).V(2)->SetS();
			*Nm=*Nm+(*fi2).N();
			k++;
		}
	}
	*Nm=*Nm/k;
	node->m=m;
}

inline minCFT BestAttachPosition(ToothNode *node,Point3f F,Point3f T ,double base[3][3],bool pickarea,int youxianji = 0)//面片上最佳位置的选取//20170906修改
{ 
	// Firstly,choose the AttachPosition by 
	//m is the model of tooth,including toothroot and crown
	//attach is the model of attachment
	//attachface is the model of attachmentface
	//facePoint is a matrix whose each single column is Point index

	//float Ymin=100;
	//float Ymax=-100;
	//float Zmax=-100;
	//Point3f Csi=Cs(node->m);
	Point3f Csi=Cs(node);//新的CR点的计算函数
	node->CR=Csi;
	////find the boundary of tooth,Xb=0,Yboundmin,Yboundmax,Zboundmin,Zboundmax
	//CMeshO::VertexIterator Vi;
	//int K=0;
	//for (Vi = m->cm.vert.begin();Vi!=m->cm.vert.end();Vi++)
	//{  
	//	K++;
	//	CVertexO *V = &(*Vi);
	//	vcg::Point3f vertex=V->P();
	//	vertex=M(base,vertex-Csi);//进行坐标变换
	//	if (vertex.Y()>Ymax)//进行边界的选择
	//	{
	//		Ymax=vertex.Y();
	//	}
	//	if (vertex.Y()<Ymin)
	//	{
	//		Ymin=vertex.Y();
	//	}
	//	if (vertex.Z()>Zmax)
	//	{
	//		Zmax=vertex.Z();
	//	}
	//} 
	//float Xbound=0;
	//float Yboundmin=2*Ymin/5.0;//Yboundmin
	//float Yboundmax=2*Ymax/5.0;//Yboundmax
	//float Zboundmin=Zmax/5.5;//Zboundmin
	//float Zboundmax=2.5*Zmax/4.0;//Zboundmax
	//cout<<Ymin<<"    "<<Ymax<<"    "<<Zmin<<"    "<<Zmax<<"    "<<endl;//test the boundaries
	Point3f Nm(0,0,0);
	Point3f TR=M(node->toothFeature.initlateralPivot,node->barycentric_coord,node->pca[1],node->pca[2],node->pca[0]);
	float Ymin=-0.5*node->toothFeature.mesialDistalWidth;
	//近舌点作为Ymax
	float Ymax=0.5*node->toothFeature.mesialDistalWidth;

	float Zmax=(M(node->toothFeature.initlateralCusp,node->barycentric_coord,node->pca[1],node->pca[2],node->pca[0])).Z();
	float Zmin=(M(node->toothFeature.initlateralGump,node->barycentric_coord,node->pca[2],node->pca[2],node->pca[0])).Z();

	float Zboundmin=Zmin/3;//Zboundmin
	float Xbound=0;
	float Yboundmin=2.5*Ymin/5.0;//Yboundmin
	float Yboundmax=2.5*Ymax/5.0;//Yboundmax
	float Zboundmax=2.0*Zmax/4.0;//zboundmax
	ChooseAttachFace(node,&Nm,Xbound,Yboundmin,Yboundmax,Zboundmin,Zboundmax,base,TR,pickarea);
	MeshModel *attachface=node->m;//we start copy the attachface from m,in order to delete these unsatisfied faces
	//float kr=1.0;
	CMeshO::FaceIterator fi;
	/*	CMeshO::VertexIterator Vi;*/
	CMeshO::VertexIterator V2i;
	int j=0;
	float Tmin=1000000000;
	int minind;
	Point3f minCoor;
	int k=0;
	Point3f Fac;
	Point3f d;
	minCFT final;
	tri::UpdateNormals<CMeshO>::PerVertex(attachface->cm);
	for (V2i = attachface->cm.vert.begin();(V2i!=attachface->cm.vert.end());V2i++)//face.size
	{   
		if((*V2i).IsS() ){
			k++;
			CVertexO *V2 = &(*V2i);
			vcg::Point3f d=V2->P();
			d-=Csi;
			d=d/d.Norm()*(d.Norm()+node->doa);
			Point3f Fac = F;

			float mis;
			if(youxianji == 4)//其他变换
				Fac=OptimFT(d,T,F);

			else if(youxianji == 0) //绕Z轴旋转
			{
				Fac=OptimFT(d,T,F,-Nm,pickarea);//与Nm垂直

			}
			else if(youxianji == 1){//绕X轴旋转
				Fac=OptimFT(d,T,F,Nm);
				//mis=(Fac-F).SquaredNorm()+((d^Fac)-T).SquaredNorm()+lamda*Fac.SquaredNorm();
			}else if(youxianji == 2){//绕y轴旋转
				Fac=OptimFT(d,T,F,-Nm,pickarea);
			}
			else if(youxianji == 3){//平移
				Fac=OptimFT(d,T,F);
			}
			mis=(Fac-F).SquaredNorm()+((d^Fac)-T).SquaredNorm()+lamda*Fac.SquaredNorm();
			//cout << "Mis:" << mis << endl;
			if (mis<Tmin)
			{
				Tmin=mis;//相对误差
				minind=j;
				minCoor=d+Csi;
				final.Fa=Fac;
				final.Ta=d^Fac;
				final.bp=minCoor;
			}
		}
		j++;
	}
	return final;
};


inline minCFT BestAttachPosition(MeshModel *m,Point3f F,Point3f T ,double base[3][3],bool pickarea,int youxianji = 0)//面片上最佳位置的选取
{ 
	MeshModel *attachface=m;
	float Ymin=100;
	float Ymax=-100;
	float Zmax=-100;
	Point3f Csi=Cs(m);
	//find the boundary of tooth,Xb=0,Yboundmin,Yboundmax,Zboundmin,Zboundmax
	CMeshO::VertexIterator Vi;
	int K=0;
	for (Vi = m->cm.vert.begin();Vi!=m->cm.vert.end();Vi++)
	{  
		K++;
		CVertexO *V = &(*Vi);
		vcg::Point3f vertex=V->P();
		vertex=M(base,vertex-Csi);//进行坐标变换
		if (vertex.Y()>Ymax)//进行边界的选择
		{
			Ymax=vertex.Y();
		}
		if (vertex.Y()<Ymin)
		{
			Ymin=vertex.Y();
		}
		if (vertex.Z()>Zmax)
		{
			Zmax=vertex.Z();
		}
	} 
	float Xbound=0;
	float Yboundmin=2*Ymin/5.0;//Yboundmin
	float Yboundmax=2*Ymax/5.0;//Yboundmax
	float Zboundmin=Zmax/5.5;//Zboundmin
	float Zboundmax=2.5*Zmax/4.0;//Zboundmax
	Point3f Nm(0,0,0);

	ChooseAttachFace(attachface,&Nm,Xbound,Yboundmin,Yboundmax,Zboundmin,Zboundmax,base,Csi,pickarea);
	CMeshO::FaceIterator fi;
	CMeshO::VertexIterator V2i;
	int j=0;
	float Tmin=1000000000;
	int minind;
	Point3f minCoor;
	int k=0;
	Point3f Fac;
	Point3f d;
	minCFT final;
	tri::UpdateNormals<CMeshO>::PerVertex(attachface->cm);
	for (V2i = attachface->cm.vert.begin();(V2i!=attachface->cm.vert.end());V2i++)//face.size
	{   
		if((*V2i).IsS() ){
			k++;
			CVertexO *V2 = &(*V2i);
			vcg::Point3f d=V2->P();
			d-=Csi;
			Point3f Fac = F;

			float mis;
			if(youxianji == 4)//其他变换
				Fac=OptimFT(d,T,F);

			else if(youxianji == 0) //绕Z轴旋转
			{
				Fac=OptimFT(d,T,F,-Nm,pickarea);
			}
			else if(youxianji == 1){//绕X轴旋转
				Fac=OptimFT(d,T,F,Nm);
			}else if(youxianji == 2){//绕y轴旋转
				Fac=OptimFT(d,T,F,-Nm,pickarea);
			}
			else if(youxianji == 3){//平移
				Fac=OptimFT(d,T,F);
			}
			mis=(Fac-F).SquaredNorm()+((d^Fac)-T).SquaredNorm()+lamda*Fac.SquaredNorm();
			if (mis<Tmin)
			{
				Tmin=mis;//相对误差
				minind=j;
				minCoor=d+Csi;
				final.Fa=Fac;
				final.Ta=d^Fac;
				final.bp=minCoor;
			}
		}
		j++;
	}
	return final;
};

void ForceEditPlugin::Calattloc(bool anchorage,Point3f Ft)
{//计算牙齿附件的函数模块//20170906修改
	//修改
    ofstream fp("errorof3e.txt");
	foreach(ToothNode *tni, teeth_tree[current_teeth]->seqtooth)
	{//使用自定义的foreach,对当前牙齿的seqtooth进行遍历,初始化使用tni
		std::string s=tni->nodename;
		tni->attforces.clear();
		tni->attloc.clear();
		double B[3][3]=
		{
			{tni->init_pca[1].X(),tni->init_pca[2].X(),tni->init_pca[0].X()},
			{tni->init_pca[1].Y(),tni->init_pca[2].Y(),tni->init_pca[0].Y()},
			{tni->init_pca[1].Z(),tni->init_pca[2].Z(),tni->init_pca[0].Z()}
		};
		double base[3][3];
		GetMatrixInverse(3,B,base);//逆矩阵的计算是正确的，计算得到逆矩阵base
		minCFT final, final2;
		//cout<<"final"<<"*******************"<<endl;
		if(tni->mainloc >=0)
		{
			cout << tni->nodename << endl;
			if(tni->mainloc >=3){///牙齿做平移移动
				final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,3);//加了
			/*	cout<<"final被赋值"<<endl;*/
				tni->cft=final;
				tni->attforces.push_back(tni->n_TargetForce);
				tni->attloc.push_back(final.bp);
				/*	cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;*/
			}
			else if(tni->mainloc == 0){//绕Z轴旋转
				final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,0);//计算第一次迭代所产生的力系
				cout<<"final被赋值"<<endl;
				tni->cft=final;
				tni->attforces.push_back(final.Fa);//将力和力矩存入tni这个类中
				tni->attloc.push_back(final.bp);
				//Point3f F2=tni->n_TargetForce-final.Fa;//减去第一个位置优化产生的力系
				//Point3f T2=tni->n_TargetTorque-final.Ta;
				//minCFT final2 = BestAttachPosition(tni,F2,T2,base,false,0);//计算第二次迭代所产生的力系
				//            tni->cft2=final2;
				//tni->attforces.push_back(final2.Fa);//力系的存储
				//tni->attloc.push_back(final2.bp);
	/*			0"牙齿编号："<<s<<"旋转类型:绕z旋转"<<endl;
				cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;*/
				/*	cout<<"Force2  "<<final2.Fa.X()<<"  "<<final2.Fa.Y()<<"  "<<final2.Fa.Z()<<endl;*/
			}
			else if(tni->mainloc == 1)
			{///绕X轴旋转

				if (s =="UL4"||s=="UL5"||s=="UR4"||s=="UR5"||s=="LL4"||s=="LL5"||s=="LR4"||s=="LR5")
				{
					final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,1);//加了
                    cout<<"final被赋值"<<endl;
					tni->attforces.push_back(final.Fa);
					tni->attloc.push_back(final.bp);
					Point3f F2=tni->n_TargetForce-final.Fa;//减去第一个位置优化产生的力系
					Point3f T2=tni->n_TargetTorque-final.Ta;
					//minCFT finla2 = BestAttachPosition(tni->m,-tni->attforces[0],tni->n_TargetTorque,base,false,1);
					final2 = BestAttachPosition(tni,F2,T2,base,true,1);
					tni->attforces.push_back(final2.Fa);
					tni->attloc.push_back(final2.bp);
	/*				cout<<"牙齿编号："<<s<<"旋转类型:绕x旋转"<<endl;
					cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;
					cout<<"Force2  "<<final2.Fa.X()<<"  "<<final2.Fa.Y()<<"  "<<final2.Fa.Z()<<endl;*/
				}else{
					final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,1);//加了
                    
					tni->attforces.push_back(final.Fa);
					tni->attloc.push_back(final.bp);
					//Point3f F2=tni->n_TargetForce-final.Fa;//减去第一个位置优化产生的力系
					//Point3f T2=tni->n_TargetTorque-final.Ta;
					////minCFT finla2 = BestAttachPosition(tni->m,-tni->attforces[0],tni->n_TargetTorque,base,false,1);
					//minCFT final2 = BestAttachPosition(tni,F2,T2,base,true,1);
					//tni->attforces.push_back(final2.Fa);
					//tni->attloc.push_back(final2.bp);
					//cout<<"牙齿编号："<<s<<"旋转类型:绕x旋转"<<endl;
					//cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;
					//cout<<"Force2  "<<final2.Fa.X()<<"  "<<final2.Fa.Y()<<"  "<<final2.Fa.Z()<<endl;
				}

			}	
			else if(tni->mainloc == 2)
			{///绕Y轴旋转
				//bool q=s=="UL1";
				if (s=="LL1"||s=="LR1"||s=="LL2"||s=="LR2"||s=="UL1"||s=="UR1"||s=="UL2"||s=="UR2")
				{
					final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,2);//加了

					tni->attforces.push_back(final.Fa);
					tni->attloc.push_back(final.bp);
					Point3f F2=tni->n_TargetForce-final.Fa;
					Point3f T2=tni->n_TargetTorque-final.Ta;
					final2 = BestAttachPosition(tni,F2,T2,base,false,2);

					tni->attforces.push_back(final2.Fa);
					tni->attloc.push_back(final2.bp);
		/*			cout<<"牙齿编号："<<s<<"旋转类型:绕y旋转"<<endl;
					cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;
					cout<<"Force2  "<<final2.Fa.X()<<"  "<<final2.Fa.Y()<<"  "<<final2.Fa.Z()<<endl;*/
				}else{
					final=BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,true,2);//加了

					tni->attforces.push_back(final.Fa);
					tni->attloc.push_back(final.bp);
					//Point3f F2=tni->n_TargetForce-final.Fa;
					//Point3f T2=tni->n_TargetTorque-final.Ta;
					//minCFT final2 = BestAttachPosition(tni,F2,T2,base,false,2);

					//tni->attforces.push_back(final2.Fa);
					//tni->attloc.push_back(final2.bp);
					//cout<<"牙齿编号："<<s<<"旋转类型:绕y旋转"<<endl;
					//cout<<"Force1  "<<final.Fa.X()<<"  "<<final.Fa.Y()<<"  "<<final.Fa.Z()<<endl;
					//cout<<"Force2  "<<final2.Fa.X()<<"  "<<final2.Fa.Y()<<"  "<<final2.Fa.Z()<<endl;
				}
			}	
		}
		if(anchorage&&(tni->mainloc <0)&&((tni->id)%7==0||(tni->id)%7==6))
		{//满足条件的支抗附件的计算
			//Point3f T1,F2;
			//T1.X()=0;
			//T1.Y()=0;
			//T1.Z()=0;
			//F2.X()=Ft.X();
			//F2.Y()=Ft.Y();
			//F2.Z()=Ft.Z();
			//minCFT final = BestAttachPosition(tni,F2,T1,base,true,3);//加了
			final = BestAttachPosition(tni,tni->n_TargetForce,tni->n_TargetTorque,base,0,3);//取背面的点
			tni->attforces.push_back(final.Fa);
			tni->attloc.push_back(final.bp);
		/*	cout<<"牙齿编号"<<tni->nodename<<"支抗附件的位置"<<final.bp.X()<<"  "<<final.bp.Y()<<"  "<<final.bp.Z()<<endl;*/
		}
		Point3f n_TF_Normal = tni->n_TargetForce.Normalize();
		Point3f n_TT_Normal = tni->n_TargetTorque.Normalize();
		Point3f AAF_Normal = (final.Fa + final2.Fa).Normalize();
		Point3f AAT_Normal = (final.Ta + final2.Ta).Normalize();
		float error1 = (tni->n_TargetForce).Norm() == 0 ? 0 : acos(n_TF_Normal * AAF_Normal);
		float error2 = (tni->n_TargetTorque).Norm() == 0 ? 0 : acos(n_TT_Normal * AAT_Normal);
	/*	cout<<"error1"<<error1<<"  "<<"error2"<<error2<<endl;*/
		fp<<tni->nodename<<" "<<error1<<" "<<error2<<endl;
	}

	hasnewfata = true;//绘制力和圆锥的flag，如果true则绘制
    fp.close();
}

///计算附件的施力的主方向
///1.当牙齿绕Z轴的旋转量超过一定程度时，需要添加附件作用(Rz)
///2.当牙齿进行近远中移动时，需要添加附件(Ty)
///3.当牙齿进行伸长压低时，需要添加附件(Tz)
///
void ForceEditPlugin::slotCalDesiredForces(){//跟踪想要的牙齿的力
	//修改

	//const float rothreh = 20.0;
	//const float tranthreh = 2.0;
	const float rothreh = 20.0;


	//定义优先级的代码
	Point3f FB;//牙套所给与的力
	FB.X()=0;
	FB.Y()=0;
	FB.Z()=0;
	float Ft;//支抗牙所给的力的数值大小
	int ina = -2;//记录附件的数目20180626
	foreach(ToothNode *tni, teeth_tree[current_teeth]->seqtooth)
		//这个不能理解啊，为什么初始化tni里面就有值了呢？foreach是个重新定义的函数，将seqtooth定义到tni里面了
	{//对teeth_tree中的牙齿进行循环

		float tranthrehx=2,tranthrehy=2,tranthrehz=2,tranthreh = 2;
		Matrix44f maintr;
		maintr.SetIdentity();
		tni->mainloc = -1;//初始化为-1
		//计算当前步骤的牙齿移动的力和力矩
		tni->C_TargetForce=Force(tni->m,tni->move_state_matrix[teeth_tree[current_teeth]->current_step+1])-Force(tni->m,tni->move_state_matrix[teeth_tree[current_teeth]->current_step]);//当前的步骤0+1，日后可以换
		/*	cout<<"单步的力"<<tni->C_TargetForce.X()<<" "<<tni->C_TargetForce.Y()<<" "<<tni->C_TargetForce.Z()<<endl;*/
		tni->C_TargetTorque=Torque(tni,tni->move_state_matrix[teeth_tree[current_teeth]->current_step+1])-Torque(tni,tni->move_state_matrix[teeth_tree[current_teeth]->current_step]);
		string s=tni->nodename;
		if (s=="UL7"||s=="UR7"||s=="UL6"||s=="UR6"||s=="LL7"||s=="LR7"||s=="LL6"||s=="LR6")
		{
			tranthrehy=0.5;

		}
		if(math::Abs(tni->rotation_theta) > rothreh){//z旋转的角度大于限定角度
			tni->mainloc = 0;
			maintr = CalRotation(tni,0,tni->init_barycentric_coord,tni->rotation_theta);
		}
		if(math::Abs(tni->x_rotataion) > rothreh){//x旋转的角度大于限定角度
			tni->mainloc = 1;
			maintr = CalRotation(tni,1,tni->init_barycentric_coord,tni->x_rotataion);
		}
		if(math::Abs(tni->y_rotataion) > rothreh){//y旋转的角度大于限定角度
			tni->mainloc = 2;
			maintr = CalRotation(tni,2,tni->init_barycentric_coord,tni->y_rotataion);
		}
		if(math::Abs(tni->x_tran) > tranthrehx){//x平移的距离大于限定的距离
			tni->mainloc = 3;
			maintr = CalTran(tni,1,tni->x_tran);
		}
		if(math::Abs(tni->y_tran) > tranthrehy){//y平移的距离大于限定的距离
			tni->mainloc = 4;
			maintr = CalTran(tni,2,tni->y_tran);
		}
		if(math::Abs(tni->extint) > tranthrehz){//z平移的距离大于限定的距离
			tni->mainloc = 5;
			maintr = CalTran(tni,0,tni->extint);
		}

		if(tni->mainloc == -1){//没有找到对应的牙齿或者运动的方式
			std::cout << "No attaachment" << std::endl;
			Matrix44f maintr1,maintr2;
			maintr1.SetIdentity();
			maintr2.SetIdentity();
			maintr1 = CalTran(tni,tni->x_tran,tni->y_tran,tni->extint);//计算平移的旋转矩阵，严格遵守x,y,z的输入
			maintr2=CalRotation(tni,tni->init_barycentric_coord,tni->x_rotataion,tni->y_rotataion,tni->rotation_theta);
			tni->n_TargetForce = Force(tni->m,maintr1);//计算平移的力
			tni->n_TargetTorque = Torque(tni,maintr2);//计算旋转的力矩

		}else if(tni->mainloc >=0 &&tni->mainloc <3){//旋转的角度很大的情况下
			std::cout << "Need double attachments " << std::endl;
			tni->n_TargetForce.SetZero();//忽略牙齿的作用力
			//TargetToque = Torque(CMN->m,tranre *rf * trano);
			tni->n_TargetTorque = Torque(tni,maintr);//忽略牙齿的作用力
			//tni->n_TargetTorque.Normalize();//标准化
			ina++;
		}else{
			std::cout << "Need single attaachment" << std::endl;
			tni->n_TargetForce = Force(tni->m,maintr);//平移的距离很大的话
			if (tni->mainloc==4 && (s=="UL7"||s=="UL6"||s=="UL5"||s=="UL4"||s=="UL3"||s=="UL2"||s=="UL1"))//局部坐标系的y轴反了
			{
				tni->n_TargetForce = -tni->n_TargetForce;
			}
			//tni->n_TargetForce.Normalize();//标准化 
			//TargetToque = Torque(CMN->m,tranre *rf * trano);
			tni->n_TargetTorque.SetZero();
			ina++;
		}

		Point3f TF = tni->n_TargetForce;//计算每一颗牙齿的力
		/*	cout<<"编号是"<<tni->nodename<<"力的大小是"<<TF.Norm()<<"分别是"<<TF.X()<<","<<TF.Y()<<","<<TF.Z()<<endl;*/

		if (!((tni->id)%7==0||(tni->id)%7==6))
		{
			FB+=TF;
		}

	}
	FB/=4;//总共有4颗牙齿作为支抗牙，要求解得到每颗牙齿的大小
	float maxFt=100.00;
	Ft = FB.Norm();
	if (Ft < maxFt || ina >= 6)//如果小于移动的力的阈值或者有附件的牙齿数目大于6
	{
		teeth_tree[current_teeth]->hasanchorage=0;
	}else{
		//如果达到移动的力的阈值
		teeth_tree[current_teeth]->hasanchorage=1;
	}
	cout<<"是否需要支抗附件:"<<teeth_tree[current_teeth]->hasanchorage<<endl;
	cout << "牙齿" <<"  "<<ina << endl;
	Calattloc(teeth_tree[current_teeth]->hasanchorage,FB);//计算附件的位置
	forceJudge=true;
}


void ForceEditPlugin::slotCalRefTorque()
{
	const float rothreh = 20.0;
	const float tranthreh = 2.0;

	vcg::Point3f FB = vcg::Point3f(0,0,0);
	float Ft;
	foreach(ToothNode* tni,teeth_tree[current_teeth]->seqtooth)
	{
		Matrix44f maintr;
		maintr.SetIdentity();
		tni->mainloc = -1;

		if(math::Abs(tni->rotation_theta) > rothreh)//z旋转的角度大于限定角度
		{
			tni->mainloc = 0;
			maintr = CalRotation(tni,0,tni->init_barycentric_coord,tni->rotation_theta);
		}
		if(math::Abs(tni->x_rotataion) > rothreh)//x旋转的角度大于限定角度
		{
			tni->mainloc = 1;
			maintr = CalRotation(tni,1,tni->init_barycentric_coord,tni->x_rotataion);
		}
		if(math::Abs(tni->y_rotataion) > rothreh)//y旋转的角度大于限定角度
		{
			tni->mainloc = 2;
			maintr = CalRotation(tni,2,tni->init_barycentric_coord,tni->y_rotataion);
		}
		if(math::Abs(tni->x_tran) > tranthreh)//x平移的距离大于限定的距离
		{
			tni->mainloc = 3;
			maintr = CalTran(tni,1,tni->x_tran);
		}
		if(math::Abs(tni->y_tran) > tranthreh)//y平移的距离大于限定的距离
		{
			tni->mainloc = 4;
			maintr = CalTran(tni,2,tni->y_tran);
		}
		if(math::Abs(tni->extint) > tranthreh)//z平移的距离大于限定的距离
		{
			tni->mainloc = 5;
			maintr = CalTran(tni,0,tni->extint);
		}

		if(tni->mainloc == -1)
		{
			qDebug()<<__FUNCTION__<<__LINE__ << " No attaachment";
			Matrix44f maintr1,maintr2;
			maintr1.SetIdentity();
			maintr2.SetIdentity();
			maintr1 = CalTran(tni,tni->x_tran,tni->y_tran,tni->extint);//计算平移的旋转矩阵，严格遵守x,y,z的输入
			maintr2=CalRotation(tni,tni->init_barycentric_coord,tni->x_rotataion,tni->y_rotataion,tni->rotation_theta);
			tni->n_TargetForce = Force(tni,maintr1);//计算平移的力
			tni->n_TargetTorque = Torque(tni,maintr2);//计算旋转的力矩
		}
		else if(tni->mainloc >=0 &&tni->mainloc <3)//旋转的角度很大的情况下
		{
			qDebug()<<__FUNCTION__<<__LINE__ << " Need double attachments";
			tni->n_TargetForce.SetZero();//忽略牙齿的作用力
			tni->n_TargetTorque = Torque(tni,maintr);//忽略牙齿的作用力
		}
		else
		{
			qDebug()<<__FUNCTION__<<__LINE__ << " Need single attaachment";
			tni->n_TargetForce = Force(tni,maintr);//平移的距离很大的话
			tni->n_TargetTorque.SetZero();
		}

		Point3f TF = tni->n_TargetForce;//计算每一颗牙齿的力
		std::cout<<__FUNCTION__<<__LINE__<<" "<<tni->nodename<<" 力的大小是 "<<TF.Norm()<<" 分别是 "<<TF.X()<<","<<TF.Y()<<","<<TF.Z()<<std::endl;

		if (!((tni->id)%7==0||(tni->id)%7==6))
		{
			FB+=TF;
		}
	}
	FB/=4;//总共有4颗牙齿作为支抗牙，要求解得到每颗牙齿的大小
	float maxFt=100.00;
	Ft=FB.Norm();
	if (Ft<maxFt)//如果小于移动的力的阈值
	{
		teeth_tree[current_teeth]->hasanchorage=0;
	}else{
		//如果达到移动的力的阈值
		teeth_tree[current_teeth]->hasanchorage=1;
	}
	std::cout<<__FUNCTION__<<__LINE__<<" 是否需要支抗附件:"<<teeth_tree[current_teeth]->hasanchorage<<std::endl;;
	Calattloc(teeth_tree[current_teeth]->hasanchorage,FB);//计算附件的位置
}

void ForceEditPlugin::slotMeasure()
{
	if(is_measure)
	{
		gla->setCursor(QCursor(QPixmap(":/images/cur_info.png"),1,1));
		is_measure = false;
	}
	else{
		was_ready = true;
		rubberband.Reset();
		is_measure = true;
		gla->setCursor(QCursor(QPixmap(":/images/cur_measure.png"),15,15));
	}
	firstclick = true;
}

//--------------------------------------------------------------------------------------------------------
void ForceEditPlugin::importAttFile()
{
	//QString dir_f = "C:\\Users\\User\\Desktop\\out.txt";
	//QByteArray ba = dir_f.toLatin1();    
	//ofstream exput;
	//exput.open(ba.data());
	//if(!exput.is_open())
	//{
	//	qDebug()<<__FUNCTION__<<__LINE__;
	//}

	QString dir = QFileDialog::getOpenFileName(this->gla, tr("Open Att File..."),md->mm()->fullName(),tr("Setting (*.att)"));
	FILE *fp = fopen(dir.toLocal8Bit().data(), "rb");
	if (fp == NULL) {  
		perror("Open file error");  
		return;  
	}

	int emrID = -1;
	int attNum = 0; // 附件个数
	fread(&emrID, sizeof(int), 1, fp);	// 病例ID
	fread(&attNum, sizeof(int), 1, fp);	// 附件个数
	printf("[%s:%d] ID: %d 附件个数: %d\n", __FUNCTION__, __LINE__, emrID, attNum);
	QString attName, parentName;
	for (int i = 0; i < attNum; ++i)
	{
		char nameBuf[50];
		int nameLen = 0;
		int attType = 0;// 附件类型
		vcg::Point2i duration; // 起始和终止阶段
		fread(&attType, sizeof(int), 1, fp);// 附件类型
		fread(&nameLen, sizeof(int), 1, fp);// 附件名称长度
		fread(nameBuf, 1, nameLen, fp);		// 附件名称
		nameBuf[nameLen] = 0;
		attName = QString::fromLocal8Bit(nameBuf);
		//exput<<attName.toLocal8Bit().data()<<"\n";
		fread(&nameLen, sizeof(int), 1, fp);// 父节点名称长度
		fread(nameBuf, 1, nameLen, fp);		// 父节点名称
		nameBuf[nameLen] = 0;
		parentName = QString::fromLocal8Bit(nameBuf);
		fread(&duration, sizeof(int), 2, fp);	// 出现的阶段a与消失的阶段b, 在a之前粘贴, 在b之后去除, 全为0表示从始至终
		//exput<<duration[0]<<" "<<duration[1]<<"\n";
		{
			int vertexNum;
			int faceNum;
			fread(&vertexNum, sizeof(int), 1, fp);	// 顶点数
			fread(&faceNum, sizeof(int), 1, fp);	// 面片数
			//exput<<vertexNum<<"\n";
			//exput<<faceNum<<"\n";
			vcg::Point3f *vertices = new vcg::Point3f[vertexNum];
			fread((char*)(vertices), sizeof(float), 3 * vertexNum, fp);
			//for (int i=0;i<vertexNum;i++)
			//{
			//	exput<<vertices[i].X()<<" "<<vertices[i].Y()<<" "<<vertices[i].Z()<<"\n";
			//}
			delete[]vertices;
			vcg::Point3i *vhandle = new vcg::Point3i[faceNum];
			fread((char*)(vhandle), sizeof(int), 3 * faceNum, fp);
			//for (int i=0; i<faceNum;i++)
			//{
			//	exput<<vhandle[i].X()<<" "<<vhandle[i].Y()<<" "<<vhandle[i].Z()<<"\n";
			//}
			//exput<<std::endl;
			delete[]vhandle;
			qDebug()<<__FUNCTION__<<__LINE__<<" "<<attType<<" "<<attName<<" "<<duration[0]<<" "<<duration[1]<<" "<<vertexNum<<" "<<faceNum;
		}
		
		AttachmentNode *tni = teeth_tree[0]->GetAttByName(attName);
		if(!tni) {
			tni = teeth_tree[1]->GetAttByName(attName);
		}
		if(!tni) {
			std::cout << "Not Find Node!!!"<< std::endl;
			continue;
		}
		ToothNode* pNi = teeth_tree[0]->GetToothByQName(parentName);
		if (!pNi)
		{
			pNi = teeth_tree[1]->GetToothByQName(parentName);
		}
		if(!pNi) {
			std::cout << "Not Find Node!!!"<< std::endl;
			continue;
		}
		
		tni->attType = attType;
		tni->duration[0] = duration.X();
		tni->duration[1] = duration.Y();
		tni->parent = pNi;
	}
	fclose(fp);
	//exput.close();

	printf("[%s:%d] Attachment file loaded [%s].\n", __FUNCTION__, __LINE__, dir.toLocal8Bit().data());
	return;
}


void ForceEditPlugin::slotImportAttInfo()
{
	QString dir = QFileDialog::getOpenFileName(this->gla, tr("Open FA File..."),md->mm()->fullName(),tr("Setting (*.fa)"));
	char*  ch;
	QByteArray ba = dir.toLatin1();
	ch=ba.data();
	ifstream inputs;
	std::string oneline;
	inputs.open(ch);
	if(inputs.is_open())
	{
		int attnum;
		getline(inputs, oneline);
		istringstream strStream(oneline);
		strStream >> attnum;
		strStream.clear();
		for(int j = 0; j< attnum ;j++)
		{
			std::string attname;
			getline(inputs, oneline);
			strStream.str(oneline);
			strStream >> attname;
			strStream.clear();
			AttachmentNode *tni = teeth_tree[0]->GetAttByName(attname);
			if(!tni) {
				tni = teeth_tree[1]->GetAttByName(attname);
			}
			if(!tni) {
				std::cout << "Not Find Node!!!"<< std::endl;
				continue;
			}
			getline(inputs, oneline);
			strStream.str(oneline);
			strStream >> tni->barycentric_coord.X() >> tni->barycentric_coord.Y() >> tni->barycentric_coord.Z();
			strStream >> tni->pca[1].X() >> tni->pca[1].Y() >> tni->pca[1].Z();
			strStream >> tni->pca[2].X() >> tni->pca[2].Y() >> tni->pca[2].Z();
			strStream >> tni->pca[0].X() >> tni->pca[0].Y() >> tni->pca[0].Z();
			strStream >> tni->mp.X() >> tni->mp.Y() >> tni->mp.Z();
			strStream >> tni->nv.X() >> tni->nv.Y() >> tni->nv.Z();
			strStream.clear();
			for(int i =0;i<3;i++){
				tni->old_pca[i] = tni->pca[i];
			}
			tni->orign_barycentric_coord = tni->barycentric_coord;
		}
	}
	inputs.close();

	importAttFile();
}

inline static int GetIndexVertex(MeshModel *m, CVertexO *p)
{
	return p-&*(m->cm.vert.begin());
}

static int ExportAttFile(std::vector<TeethTree *> meshtrees, QString dirname)
{
	int nodesize =0;
	for(int i =0;i<2;i++)
	{
		foreach(AttachmentNode *mni, meshtrees[i]->attNodes)
		{
			nodesize ++;
			QString attname = dirname + "\\" + mni->qnodename +".obj";
			int result = tri::io::ExporterOBJ<CMeshO>::Save(mni->m->cm,attname.toLocal8Bit().data(),true);
			if (!result)
			{
				qDebug()<<__FUNCTION__<<__LINE__<< " "<<"Save "<< attname<<" Successed.";
			}
		}

	}

	int p_id = 8;
	QString attbinname = dirname+"\\"+"attachment.att";
	ofstream rs(attbinname.toLocal8Bit().data(),ios::binary);
	rs.write((char*)&p_id,sizeof(int));
	rs.write((char*)&(nodesize),sizeof(int));
	for(int i =0;i<2;i++)
	{
		foreach(AttachmentNode *mni, meshtrees[i]->attNodes)
		{
			int type = mni->attType;
			rs.write((char*)&type,sizeof(int));
			int Appearindex, Disappearindex;
			Appearindex = mni->duration[0];
			Disappearindex = mni->duration[1];
			int vertexsnum,facesnum;
			facesnum = 0;
			vertexsnum = 0;
			MeshModel *mmi = mni->m;
			CMeshO::VertexIterator vi;
			for (vi = mni->m->cm.vert.begin();vi!=mni->m->cm.vert.end();vi++)
			{
				if(!(*vi).IsD())
					vertexsnum++;
			}
			CMeshO::FaceIterator fi;
			for(fi = mni->m->cm.face.begin();fi != mni->m->cm.face.end();fi++)
			{
				if(!(*fi).IsD())
					facesnum ++;
			}
			int nodenameLen = mni->nodename.size();
			rs.write((char*)&(nodenameLen),sizeof(int));
			rs.write(mni->nodename.c_str(),nodenameLen);

			int parentlen = mni->parent->nodename.size();
			rs.write((char*)&(parentlen),sizeof(int));
			rs.write(mni->parent->nodename.c_str(),parentlen);

			rs.write((char*)&(Appearindex),sizeof(int));
			rs.write((char*)&(Disappearindex),sizeof(int));

			rs.write((char*)&(vertexsnum),sizeof(int));
			rs.write((char*)&(facesnum),sizeof(int));

			cout<< __FUNCTION__<<__LINE__<< "|4B:" << type<<endl;
			cout<<__FUNCTION__<<__LINE__ << "|4B:" << Appearindex<<endl;
			cout<<__FUNCTION__<<__LINE__ << "|4B:" << Disappearindex << endl;
			cout<<__FUNCTION__<<__LINE__ << "|4B:" << mni->nodename << endl;
			cout<<__FUNCTION__<<__LINE__ << "|4B:" << vertexsnum << endl;
			cout<<__FUNCTION__<<__LINE__ << "|4B:" << facesnum << endl;
			
			std::vector<int> VertexId(mmi->cm.vert.size());
			int numvert = 0;

			for(vi = mmi->cm.vert.begin();vi!=mmi->cm.vert.end();vi++)
			{
				if (!(*vi).IsD())
				{
					VertexId[vi-mmi->cm.vert.begin()]=numvert;
					vcg::Point3f pvi = (*vi).P();
					rs.write((char*)&pvi.X(),sizeof(float));
					rs.write((char*)&pvi.Y(),sizeof(float));
					rs.write((char*)&pvi.Z(),sizeof(float));
					numvert++;
				}
			}

			for(fi = mmi->cm.face.begin();fi != mmi->cm.face.end();fi++)
			{
				if (!(*fi).IsD())
				{
					for(int k =0;k<(*fi).VN();k++)
					{
						int vInd = VertexId[GetIndexVertex(mmi, (*fi).V(k))] + 1;
						rs.write((char*)&vInd,sizeof(int));
					}
				}
			}
		}
	}
	rs.close();
	return nodesize;
}
void ForceEditPlugin::slotExportAttInfo()
{
	QString dir = QFileDialog::getExistingDirectory(this->gla, tr("Save Atts Directory.."),"",
		QFileDialog::ShowDirsOnly| QFileDialog::DontResolveSymlinks);

	if (dir.isNull()||dir.isEmpty())
	{
		return;
	}
	
	QDir dirr(dir);
	if (!dirr.exists())
	{
		bool ok = dirr.mkdir(dir);
	}
	else
	{
		dirr.setFilter(QDir::Files);
		int fileCount = dirr.count();
		for (int i = 0; i < fileCount; i++)
			dirr.remove(dirr[i]); 
	}

	int nodesize = ExportAttFile(teeth_tree,dir);

	if(nodesize > 0)
	{
		char*  ch_fa;
		QString dir_f = dir+"\\att.fa";
		QByteArray ba = dir_f.toLatin1();    
		ch_fa=ba.data();
		ofstream exput;
		exput.open(ch_fa);
		if(exput.is_open())
		{
			exput << nodesize << "\n";
			for(int i =0;i<2;i++)
			{	
				foreach(AttachmentNode *mni , teeth_tree[i]->attNodes)
				{
					exput << mni->nodename << "\n";
					exput << mni->barycentric_coord.X() << " " << mni->barycentric_coord.Y() << " " << mni->barycentric_coord.Z() << " ";
					exput << mni->pca[1].X() << " " << mni->pca[1].Y() << " " << mni->pca[1].Z() << " ";
					exput << mni->pca[2].X() << " " << mni->pca[2].Y() << " " << mni->pca[2].Z() << " ";
					exput << mni->pca[0].X() << " " << mni->pca[0].Y() << " " << mni->pca[0].Z() << " ";
					exput << mni->mp.X() << " " << mni->mp.Y() << " " << mni->mp.Z() << " ";
					exput << mni->nv.X() << " " << mni->nv.Y() << " " << mni->nv.Z() << " "<<"\n";
				}
			}  
		}
		exput << std::endl;
		exput.close();
	}
}

inline static bool Tran4to3(vcg::Matrix44f tr,vcg::Point3f &tran,vcg::Matrix33f &res)
{
	for(int i =0;i < 3;i++){
		for(int j =0;j<3;j++){
			res[i][j] = tr[i][j];
		}
	}
	tran.X() = res[3][0];tran.Y() = res[3][1];tran.Z() = res[3][2];
	return 1;
}

static void GumDeformationByLaplacian(ToothNode *gumm, ToothNode *CMN, vcg::Matrix44f r){
	int i,j,k;
	CMeshO::VertexIterator vi;
	int degree;//顶点度数
	int v_weight = 100;//权重，越大点越固定
	int VN;//处理顶点的数量
	int flag;
	CMeshO::PerVertexAttributeHandle<int> col = tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(gumm->m->cm ,std::string("v_col"));//给每个顶点增加col属性，以便快速找到该顶点在矩阵中的列号
	for(vi=gumm->m->cm.vert.begin(),i=0; vi!=gumm->m->cm.vert.end(); ++vi,++i)
	{
		if (!(*vi).IsD())
		{				
			col[vi] = -1;//初始化列号属性
		}
	}
	VN=CMN->originhandleVex.size()+CMN->origindeformationVex.size()+CMN->originstaticVex.size();
	j = 0;
	for (i = 0;i < CMN->originhandleVex.size();i++,j++)
	{
		flag = CMN->originhandleVex[i].id;
		vi = gumm->m->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}
	for (i = 0;i < CMN->origindeformationVex.size();i++,j++)
	{
		flag = CMN->origindeformationVex[i].id;
		vi = gumm->m->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}
	for (i = 0;i < CMN->originstaticVex.size();i++,j++)
	{
		flag = CMN->originstaticVex[i].id;
		vi = gumm->m->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}

	vector<Eigen::Triplet<float>> trips;
	j = 0;
	for(i=0;i<CMN->originhandleVex.size();++i,++j)//构造矩阵
	{
		vi = gumm->m->cm.vert.begin();
		flag = CMN->originhandleVex[i].id;
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Eigen::Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Eigen::Triplet<float>(j, j,degree));
		}
	}
	for(i=0;i<CMN->origindeformationVex.size();++i,++j)//构造矩阵
	{
		vi = gumm->m->cm.vert.begin();
		flag = CMN->origindeformationVex[i].id;
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Eigen::Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Eigen::Triplet<float>(j, j,degree));
		}
	}
	for(i=0;i<CMN->originstaticVex.size();++i,++j)//构造矩阵
	{
		vi = gumm->m->cm.vert.begin();
		flag = CMN->originstaticVex[i].id;
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Eigen::Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Eigen::Triplet<float>(j, j,degree));
		}
	}
	Eigen::SparseMatrix<float> A0;//Laplacian方阵
	A0.resize(VN,VN);
	A0.setFromTriplets(trips.begin(), trips.end());

	//构造初始坐标向量
	Eigen::VectorXf X0;
	X0.resize(VN);
	Eigen::VectorXf Y0;
	Y0.resize(VN);
	Eigen::VectorXf Z0;
	Z0.resize(VN);
	j = 0;
	for(i =0;i<CMN->originhandleVex.size();i++,j++)
	{
		X0(j)=CMN->originhandleVex[i].cv.X();
		Y0(j)=CMN->originhandleVex[i].cv.Y();
		Z0(j)=CMN->originhandleVex[i].cv.Z();
	}
	for(i =0;i<CMN->origindeformationVex.size();i++,j++)
	{
		X0(j)=CMN->origindeformationVex[i].cv.X();
		Y0(j)=CMN->origindeformationVex[i].cv.Y();
		Z0(j)=CMN->origindeformationVex[i].cv.Z();
	}
	for(i =0;i<CMN->originstaticVex.size();i++,j++)
	{
		X0(j)=CMN->originstaticVex[i].cv.X();
		Y0(j)=CMN->originstaticVex[i].cv.Y();
		Z0(j)=CMN->originstaticVex[i].cv.Z();
	}
	Eigen::VectorXf BX0;
	BX0.resize(VN);
	BX0 = A0*X0;
	Eigen::VectorXf BY0;
	BY0.resize(VN);
	BY0 = A0*Y0;
	Eigen::VectorXf BZ0;
	BZ0.resize(VN);
	BZ0 = A0*Z0;

	int VN1=VN+CMN->originhandleVex.size()+CMN->originstaticVex.size();
	Eigen::VectorXf BX1;
	BX1.resize(VN1);
	Eigen::VectorXf BY1;
	BY1.resize(VN1);
	Eigen::VectorXf BZ1;
	BZ1.resize(VN1);
	for (i = 0;i<VN;i++)
	{
		BX1(i)=BX0(i);
		BY1(i)=BY0(i);
		BZ1(i)=BZ0(i);
	}
	int temp = VN;

	vcg::Matrix33f re3;
	vcg::Point3f tran;
	for(int i =0;i<3;i++)
		for(int j =0;j<3;j++)
			re3[i][j] = r[i][j];
	tran.X() = r[0][3];tran.Y() = r[1][3];tran.Z() = r[2][3];

	vi=gumm->m->cm.vert.begin();
	for(i = 0;i < CMN->originhandleVex.size();i++,temp++)
	{
		flag = CMN->originhandleVex[i].id;		
		trips.push_back(Eigen::Triplet<float>(temp, col[vi+flag], v_weight));
		Point3f p3 = (tran + re3 *(CMN->originhandleVex[i].cv))*v_weight;
		BX1(temp) = p3.X();
		BY1(temp) = p3.Y();
		BZ1(temp) = p3.Z();
	}
	for(i = 0;i < CMN->originstaticVex.size();i++,temp++)
	{
		flag = CMN->originstaticVex[i].id;		
		trips.push_back(Eigen::Triplet<float>(temp, col[vi+flag], v_weight));
		Point3f p3 = (CMN->originstaticVex[i].cv)*v_weight;
		BX1(temp) = p3.X();
		BY1(temp) = p3.Y();
		BZ1(temp) = p3.Z();
	}

	Eigen::SparseMatrix<float> A;//构造Laplacian增广矩阵
	A.resize(VN1,VN);
	A.setFromTriplets(trips.begin(), trips.end());
	Eigen::SparseMatrix<float> AT = A.transpose();
	Eigen::SparseMatrix<float> M = AT*A;

	Eigen::VectorXf BX2 = AT*BX1;
	Eigen::VectorXf BY2 = AT*BY1;
	Eigen::VectorXf BZ2 = AT*BZ1;

	Eigen::SimplicialLLT<Eigen::SparseMatrix<float>, Eigen::Lower> solverA;
	solverA.compute(M);
	Eigen::VectorXf X;
	Eigen::VectorXf Y;
	Eigen::VectorXf Z;
	X = solverA.solve(BX2);
	Y = solverA.solve(BY2);
	Z = solverA.solve(BZ2);

	j = 0;
	for (i = 0;i<CMN->originhandleVex.size();i++,j++)
	{
		CMN->handleVex[i].cv.X() = X[j];
		CMN->handleVex[i].cv.Y() = Y[j];
		CMN->handleVex[i].cv.Z() = Z[j];
	}
	for (i = 0;i<CMN->origindeformationVex.size();i++,j++)
	{
		CMN->deformationVex[i].cv.X() = X[j];
		CMN->deformationVex[i].cv.Y() = Y[j];
		CMN->deformationVex[i].cv.Z() = Z[j];
	}
	for (i = 0;i<CMN->originstaticVex.size();i++,j++)
	{
		CMN->staticVex[i].cv.X() = X[j];
		CMN->staticVex[i].cv.Y() = Y[j];
		CMN->staticVex[i].cv.Z() = Z[j];
	}
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(gumm->m->cm ,std::string("v_col"));

	vi=gumm->m->cm.vert.begin();
	vector<MarkPoint>::iterator it;
	for(it = CMN->handleVex.begin();it!= CMN->handleVex.end();++it)
	{
		int flag =it->id; 
		vi[flag].P() = it->cv;
	}
	for(it = CMN->deformationVex.begin();it!= CMN->deformationVex.end();++it)
	{
		int flag =it->id; 
		vi[flag].P() = it->cv;
	}
	for(it = CMN->staticVex.begin();it!= CMN->staticVex.end();++it)
	{
		int flag =it->id; 
		vi[flag].P() = it->cv;
	}
}

static bool ExportDevidedSTLModel(TeethTree *meshtree,const std::string pathname,bool secgum,int stepi,bool binary = true,bool savePontic = false)
{
	qDebug()<<__FUNCTION__<<__LINE__<<" "<<savePontic;
	if(secgum){
		foreach(ToothNode *mni,meshtree->nodeList){
			if(!mni->isgum) continue;
			if (stepi>meshtree->planstep-1)
			{
				RTBGum::GumDeformationByRealTimeUpdate(meshtree,mni,meshtree->initseqtooth,meshtree->planstep-1);
			}
			else
				RTBGum::GumDeformationByRealTimeUpdate(meshtree,mni,meshtree->initseqtooth,stepi);
			std::string filename = pathname + "//" + mni->nodename + ".stl";
			std::cout << __FUNCTION__<<__LINE__<<"Export " << filename << "..." << std::endl;
			FILE *fp;
			fp = fopen(filename.c_str(),"wb");
			if(fp==0)
				return 1;
			vcg::Matrix44f tr = mni->nodetr;
			vcg::Matrix33f tr3;
			vcg::Point3f tran;
			if(!Tran4to3(tr,tran,tr3))
				return 2;
			MeshModel *mm = mni->m;
			if(binary){
				// Write Header
				char header[128]="VCG                                                                                                  ";

				if(mni->nodename.length())	strncpy(header,mni->nodename.c_str(),80);
				fwrite(header,80,1,fp);
				// write number of facets
				fwrite(&(mm->cm.fn),1,sizeof(int),fp); 
				vcg::Point3f p;
				unsigned short attributes=0; //记录颜色等附加信息，暂时没有用
				CMeshO::FaceIterator fi;
				for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
				{
					// For each triangle write the normal, the three coords and a short set to zero
					p.Import(vcg::NormalizedNormal(*fi));
					p = tr3*p;;//乘以当前变换矩阵
					fwrite(p.V(),3,sizeof(float),fp);

					for(int k=0;k<3;++k){
						p=tr*(*fi).V(k)->cP();
						fwrite(p.V(),3,sizeof(float),fp);
					}
					fwrite(&attributes,1,sizeof(short),fp);
				}
			}
			else{
				if(mni->nodename.length()) fprintf(fp,"solid %s\n",mni->nodename);
				else fprintf(fp,"solid vcg\n");

				vcg::Point3f p;
				CMeshO::FaceIterator fi;	
				for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
				{
					// For each triangle write the normal, the three coords and a short set to zero
					p.Import(vcg::NormalizedNormal(*fi));
					p = tr3*p;;//乘以当前变换矩阵
					fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
					fprintf(fp,"    outer loop\n");
					for(int k=0;k<3;++k){
						p=tr*(*fi).V(k)->cP();
						fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
					}
					fprintf(fp,"    endloop\n");
					fprintf(fp,"  endfacet\n");
				}
				fprintf(fp,"endsolid vcg\n");
			}
			fclose(fp);
		}
	}
	else{
		foreach(ToothNode *mni,meshtree->nodeList){
			if(!mni->isgum) continue;
			std::string filename = pathname + "//" + mni->nodename + ".stl";
			std::cout << "Export " << filename << "..." << std::endl;
			FILE *fp;
			fp = fopen(filename.c_str(),"wb");
			if(fp==0)
				return 1;
			vcg::Matrix44f tr = mni->nodetr;
			vcg::Matrix33f tr3;
			vcg::Point3f tran;
			if(!Tran4to3(tr,tran,tr3))
				return 2;
			MeshModel *mm = mni->m;
			foreach(ToothNode *mnii,meshtree->seqtooth){
				GumDeformationByLaplacian(mni,mnii,mnii->nodetr);
			}
			//GumDeformationByRealTimeUpdate(teeths,)
			if(binary){
				// Write Header
				char header[128]="VCG                                                                                                  ";

				if(mni->nodename.length())	strncpy(header,mni->nodename.c_str(),80);
				fwrite(header,80,1,fp);
				// write number of facets
				fwrite(&(mm->cm.fn),1,sizeof(int),fp); 
				vcg::Point3f p;
				unsigned short attributes=0; //记录颜色等附加信息，暂时没有用
				CMeshO::FaceIterator fi;
				for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
				{
					// For each triangle write the normal, the three coords and a short set to zero
					p.Import(vcg::NormalizedNormal(*fi));
					p = tr3*p;;//乘以当前变换矩阵
					fwrite(p.V(),3,sizeof(float),fp);

					for(int k=0;k<3;++k){
						p=tr*(*fi).V(k)->cP();
						fwrite(p.V(),3,sizeof(float),fp);
					}
					fwrite(&attributes,1,sizeof(short),fp);
				}
			}
			else{
				if(mni->nodename.length()) fprintf(fp,"solid %s\n",mni->nodename);
				else fprintf(fp,"solid vcg\n");

				vcg::Point3f p;
				CMeshO::FaceIterator fi;	
				for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
				{
					// For each triangle write the normal, the three coords and a short set to zero
					p.Import(vcg::NormalizedNormal(*fi));
					p = tr3*p;;//乘以当前变换矩阵
					fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
					fprintf(fp,"    outer loop\n");
					for(int k=0;k<3;++k){
						p=tr*(*fi).V(k)->cP();
						fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
					}
					fprintf(fp,"    endloop\n");
					fprintf(fp,"  endfacet\n");
				}
				fprintf(fp,"endsolid vcg\n");
			}
			fclose(fp);
		} 
	}
	//输出萌出帽
	foreach(ToothNode *mni,meshtree->nodeList){
		if(mni->isgum) continue;
		if(mni->qnodename.size()<4) continue;
		if (!savePontic)
		{
			if (mni->ispontic == 1)continue;
		}
		qDebug()<<__FUNCTION__<<__LINE__<<" "<<savePontic<<" "<<mni->ispontic<<" "<<mni->isExport;
		if(savePontic && mni->ispontic == 1&&mni->isExport== false) continue;
		std::string filename = pathname + "//" + mni->nodename + ".stl";
		if (mni->iserupt)
		{
			filename = pathname + "//" + mni->nodename + ".stl";
		}
		std::cout << __FUNCTION__<<__LINE__<< "Export " << filename << "..." << std::endl;
		FILE *fp;
		fp = fopen(filename.c_str(),"wb");
		if(fp==0)
			return 1;
		vcg::Matrix44f tr = mni->nodetr;
		if (mni->showerupt)
		{
			tr = mni->eruptnode->nodetr;
		}
		vcg::Matrix33f tr3;
		vcg::Point3f tran;
		if(!Tran4to3(tr,tran,tr3))
			return 2;
		MeshModel *mm = mni->m;
		if (mni->showerupt)
		{

			mm = mni->eruptnode->m;
		}

		if(binary){
			// Write Header
			char header[128]="VCG                                                                                                  ";

			if(mni->nodename.length())	strncpy(header,mni->nodename.c_str(),80);
			fwrite(header,80,1,fp);
			// write number of facets
			fwrite(&(mm->cm.fn),1,sizeof(int),fp); 
			vcg::Point3f p;
			unsigned short attributes=0; //记录颜色等附加信息，暂时没有用
			CMeshO::FaceIterator fi;
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fwrite(p.V(),3,sizeof(float),fp);

				for(int k=0;k<3;++k){
					p=tr*(*fi).V(k)->cP();
					fwrite(p.V(),3,sizeof(float),fp);
				}
				fwrite(&attributes,1,sizeof(short),fp);
			}
		}
		else{
			if(mni->nodename.length()) fprintf(fp,"solid %s\n",mni->nodename);
			else fprintf(fp,"solid vcg\n");

			vcg::Point3f p;
			CMeshO::FaceIterator fi;	
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
				fprintf(fp,"    outer loop\n");
				for(int k=0;k<3;++k){
					p=tr*(*fi).V(k)->cP();
					fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
				}
				fprintf(fp,"    endloop\n");
				fprintf(fp,"  endfacet\n");
			}
			fprintf(fp,"endsolid vcg\n");
		}
		fclose(fp);
	}  
	//输出牙齿模型
	foreach(ToothNode *mni,meshtree->nodeList){
		if(mni->isgum) continue;
		if (mni->qnodename.size()>3) continue;
		if (!savePontic)
		{
			if (mni->ispontic == 1)continue;
		}
		qDebug()<<__FUNCTION__<<__LINE__<<" "<<savePontic<<" "<<mni->ispontic<<" "<<mni->isExport;
		if(savePontic && mni->ispontic == 1&&mni->isExport== false) continue;
		std::string filename = pathname + "//" + mni->nodename + ".stl";
		std::cout << __FUNCTION__<<__LINE__<< "Export " << filename << "..." << std::endl;
		FILE *fp;
		fp = fopen(filename.c_str(),"wb");
		if(fp==0)
			return 1;
		vcg::Matrix44f tr = mni->nodetr;
		vcg::Matrix33f tr3;
		vcg::Point3f tran;
		if(!Tran4to3(tr,tran,tr3))
			return 2;
		MeshModel *mm = mni->m;
	    if(binary){
			// Write Header
			char header[128]="VCG                                                                                                  ";

			if(mni->nodename.length())	strncpy(header,mni->nodename.c_str(),80);
			fwrite(header,80,1,fp);
			// write number of facets
			fwrite(&(mm->cm.fn),1,sizeof(int),fp); 
			vcg::Point3f p;
			unsigned short attributes=0; //记录颜色等附加信息，暂时没有用
			CMeshO::FaceIterator fi;
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fwrite(p.V(),3,sizeof(float),fp);

				for(int k=0;k<3;++k){
					p=tr*(*fi).V(k)->cP();
					fwrite(p.V(),3,sizeof(float),fp);
				}
				fwrite(&attributes,1,sizeof(short),fp);
			}
		}
		else{
			if(mni->nodename.length()) fprintf(fp,"solid %s\n",mni->nodename);
			else fprintf(fp,"solid vcg\n");

			vcg::Point3f p;
			CMeshO::FaceIterator fi;	
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
				fprintf(fp,"    outer loop\n");
				for(int k=0;k<3;++k){
					p=tr*(*fi).V(k)->cP();
					fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
				}
				fprintf(fp,"    endloop\n");
				fprintf(fp,"  endfacet\n");
			}
			fprintf(fp,"endsolid vcg\n");
		}
		fclose(fp);
	}   
	//输出附件模型
	foreach(AttachmentNode* anode,meshtree->attNodes)
	{
		//qDebug()<<__FUNCTION__<<__LINE__<<" "<<anode->duration[0]<<" "<<anode->duration[1]<<" "<<stepi;
		if(anode->duration[0] > stepi || anode->duration[1] < stepi)
		{
			qDebug()<<__FUNCTION__<<__LINE__<<" "<<anode->duration[0]<<" "<<anode->duration[1]<<" "<<stepi;
			continue;
		}
		if (( AttType::attRootType(anode->attType)==AttType::optRoot||AttType::attRootType(anode->attType)==AttType::optRot))
		{
			if (stepi == 0)//stepi == 0时输出小的
			{
				if(anode->nodename.length() == 5) continue;
			}
			else
				if (anode->nodename.length() == 4) continue;
		}
		std::string filename = pathname + "/" + anode->nodename + ".stl";
		FILE *fp;
		fp = fopen(filename.c_str(),"wb");
		if(fp==0)return 1;

		vcg::Matrix44f tr = anode->nodetr;
		vcg::Matrix33f tr3;
		vcg::Point3f tran;
		if(!Tran4to3(tr,tran,tr3)) return 2;

		MeshModel *mm = anode->m;

		if(binary)
		{
			// Write Header
			char header[128]="VCG                                                                                                  ";
			if(anode->nodename.length())	strncpy(header,anode->nodename.c_str(),80);
			fwrite(header,80,1,fp);
			// write number of facets
			fwrite(&(mm->cm.fn),1,sizeof(int),fp); 
			vcg::Point3f p;
			unsigned short attributes=0; //记录颜色等附加信息，暂时没有用
			CMeshO::FaceIterator fi;
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fwrite(p.V(),3,sizeof(float),fp);

				for(int k=0;k<3;++k){
					p=tr*(*fi).V(k)->cP();
					fwrite(p.V(),3,sizeof(float),fp);
				}
				fwrite(&attributes,1,sizeof(short),fp);
			}
		}
		else
		{
			if(anode->nodename.length()) 
				fprintf(fp,"solid %s\n",anode->nodename);
			else 
				fprintf(fp,"solid vcg\n");

			vcg::Point3f p;
			CMeshO::FaceIterator fi;	
			for(fi=mm->cm.face.begin(); fi!=mm->cm.face.end(); ++fi) if( !(*fi).IsD() )
			{
				// For each triangle write the normal, the three coords and a short set to zero
				p.Import(vcg::NormalizedNormal(*fi));
				p = tr3*p;;//乘以当前变换矩阵
				fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
				fprintf(fp,"    outer loop\n");
				for(int k=0;k<3;++k)
				{
					p=tr*(*fi).V(k)->cP();
					fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
				}
				fprintf(fp,"    endloop\n");
				fprintf(fp,"  endfacet\n");
			}
			fprintf(fp,"endsolid vcg\n");
		}
		fclose(fp);
	}

	meshtree->Gum()->m->updateDataMask(MeshModel::MM_VERTFACETOPO);
	meshtree->Gum()->m->updateDataMask(MeshModel::MM_FACEFACETOPO);

	return 0;
}

void ForceEditPlugin::slotExportModel()
{
	QString dir = QFileDialog::getExistingDirectory(this->gla, tr("Save Directory.."),"",QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
	foreach(ToothNode* nn,teeth_tree[current_teeth]->nodeList)
	{
		nn->isExport = true;
	}
	int maxStep = -1;
	qDebug()<<teeth_tree[0]->planstep << " " << teeth_tree[1]->planstep;
	switch(teeth_tree.size())
	{
	case 0:break;
	case 1:maxStep = teeth_tree[current_teeth]->planstep;break;
	case 2:maxStep = std::max(teeth_tree[0]->planstep,teeth_tree[1]->planstep);break;
	default:break;
	}
	qDebug()<<__FUNCTION__<<__LINE__<<" maxStep: "<<maxStep;
	int count = 0;
	for(int i =0;i</*teeth_tree[current_teeth]->planstep*/maxStep;i++)
	{
		count = i;
		if (i>teeth_tree[current_teeth]->planstep-1) 
			count = teeth_tree[current_teeth]->planstep-1;
		foreach(ToothNode* tnode,teeth_tree[current_teeth]->nodeList)
		{
			if (tnode->qnodename.size()>3) continue;
			if(tnode->isgum||tnode->ispontic==1){tnode->nodetr.SetIdentity();continue;}
			else if(tnode->isattachment){tnode->nodetr = tnode->parent->move_state_matrix[count];}
			else if (tnode->showerupt)
			{
				tnode->eruptnode->nodetr = tnode->eruptnode->m->cm.Tr;
				tnode->nodetr = tnode->move_state_matrix[count];
			}
			else tnode->nodetr = tnode->move_state_matrix[count];
		}
		
		foreach(AttachmentNode* anode,teeth_tree[current_teeth]->attNodes)
		{
			anode->nodetr = anode->parent->move_state_matrix[count];
		}
		
		QString fileName = dir+"/"+QString::number(i,10);
		QDir *temp = new QDir;
		if (temp->mkdir(fileName))
		{
			qDebug()<<__FUNCTION__<<__LINE__<<" "<<fileName;
			bool savePontic = forceui->ui.cb_savePontic->isChecked();
            foreach(ToothNode *tni,teeth_tree[current_teeth]->nodeList)
			{
				tni->issavem=true;
				tni->m->visible=true;
			}
			int result = ExportDevidedSTLModel(&(*teeth_tree[current_teeth]),fileName.toLocal8Bit().data(),true,i,true,savePontic);
			foreach(ToothNode *tni,teeth_tree[current_teeth]->nodeList)
			{
				tni->issavem=false;
			}
			if (result) return;
		}
		
	}
}

//--------------------------------------------------------------------------------------------------------
bool ForceEditPlugin::StartEdit(MeshDocument& _md,GLArea* _gla,TeethTree a[],bool hasprepare)
{
	
	gla = _gla;
	md = &_md;
	teeth_tree.clear();
	teeth_tree.push_back(&a[0]);//存放牙齿模型
	teeth_tree.push_back(&a[1]);//存放牙套模型
	//qDebug() << teeth_tree.size();
	gla->rm.colorMode=GLW::CMPerVert;
	gla->setCursor(QCursor(QPixmap(":/images/cur_info.png"),1,1));
	AxisSize = md->bbox().Diag() * 1.0;
	
	if(forceui == 0)
	{
		forceui = new ForceWidget(gla->window(),this);
		forceui->ui.manu_none->setChecked(true);
		connect(this,SIGNAL(suspendEditToggle()),this->gla,SLOT(suspendEditToggle()));
		connect(forceui,SIGNAL(updateMeshSetVisibilities()),this->gla,SLOT(updateMeshSetVisibilities()));
		connect(this, SIGNAL(updatemeshvisiable()), this->gla,SLOT(updateMeshSetVisibilities()));
		//其他操作
		connect(forceui->ui.pb_impAttInfo,SIGNAL(clicked()),this,SLOT(slotImportAttInfo()));//导入附件配置
		connect(forceui->ui.pb_expImqModel,SIGNAL(clicked()),this,SLOT(slotExportAttInfo()));//导出IMQ模型
 		connect(forceui->ui.pb_calTorque,SIGNAL(clicked()),this,SLOT(slotCalDesiredForces()));//计算附件参考力矩
 		connect(forceui->ui.pb_autoPaste,SIGNAL(clicked()),this,SLOT(slotAutoAtt()));//自动黏贴附件
 		connect(forceui->ui.pb_expManuModel,SIGNAL(clicked()),this,SLOT(slotExportModel()));//导出生产模型
 		connect(forceui->ui.pb_measure,SIGNAL(clicked()),this,SLOT(slotMeasure()));//测距

		//附件编辑
		connect(forceui->ui.manu_value,SIGNAL(valueChanged(double)),this,SLOT(slotGetManuValue(double)));
		connect(forceui->ui.startManu,SIGNAL(clicked()),this,SLOT(slotManStart()));
		connect(forceui->ui.applyManu,SIGNAL(clicked()),this,SLOT(slotManApply()));
		connect(forceui->ui.pb_change_D,SIGNAL(clicked()),this,SLOT(slotChangeDuration()));
		connect(forceui->ui.pb_changeType,SIGNAL(clicked()),this,SLOT(slotChangeAttType()));
		connect(forceui->ui.pushButton,SIGNAL(clicked()),this,SLOT(slotAttachOptim()));

		//附件粘贴
		QSignalMapper *attMap = new QSignalMapper;
		connect(forceui->ui.rectVertAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.rectVertAtt, AttType::rectVert);
		connect(forceui->ui.rectHoriAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.rectHoriAtt, AttType::rectHori);
		connect(forceui->ui.rectInciAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.rectInciAtt, AttType::rectInci);
		connect(forceui->ui.rectGumAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.rectGumAtt, AttType::rectGum);
		connect(forceui->ui.rectMesioAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.rectMesioAtt, AttType::rectMesio);
		connect(forceui->ui.rectDistalAtt,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.rectDistalAtt, AttType::rectDistal);
		
		connect(forceui->ui.ellipHoriAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.ellipHoriAtt, AttType::ellipHori);
		connect(forceui->ui.ellipVertAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.ellipVertAtt, AttType::ellipVert);

		connect(forceui->ui.optRotShortAtt,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.optRotShortAtt, AttType::optRotShort);
		connect(forceui->ui.optRotLongAtt,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.optRotLongAtt, AttType::optRotLong);

		connect(forceui->ui.optRootMesioAtt,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.optRootMesioAtt, AttType::optRootMesio);
		connect(forceui->ui.optRootDistalAtt,SIGNAL(clicked()),attMap,SLOT(map())); attMap->setMapping(forceui->ui.optRootDistalAtt,AttType::optRootDistal);

		connect(forceui->ui.optExtAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.optExtAtt, AttType::optExt);
		connect(forceui->ui.optDeepAtt,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.optDeepAtt, AttType::optDeep);
		//Features
		connect(forceui->ui.buttonCut_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.buttonCut_l, AttType::buttonCutout_l);
		connect(forceui->ui.buttonCut_b,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.buttonCut_b, AttType::buttonCutout);
		connect(forceui->ui.buttonCut_m,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.buttonCut_m, AttType::buttonCutout_m);
		
		connect(forceui->ui.btn_hookLeft,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.btn_hookLeft, AttType::leftHook);
		connect(forceui->ui.btn_hookLeft_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.btn_hookLeft_l, AttType::leftHook_l);
		connect(forceui->ui.btn_hookRight,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.btn_hookRight, AttType::rightHook);
		connect(forceui->ui.btn_hookRight_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.btn_hookRight_l, AttType::rightHook_l);

		connect(forceui->ui.hook_mb,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.hook_mb, AttType::hook_m_b);
		connect(forceui->ui.hook_db,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.hook_db, AttType::hook_d_b);
		connect(forceui->ui.hook_ml,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.hook_ml, AttType::hook_m_l);
		connect(forceui->ui.hook_dl,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.hook_dl, AttType::hook_d_l);

		connect(forceui->ui.powerPoint_b,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerPoint_b, AttType::powerPoint_b);
		connect(forceui->ui.powerPoint_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerPoint_l, AttType::powerPoint_l);
		connect(forceui->ui.powerRidge_b,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerRidge_b, AttType::powerRidge_b);
		connect(forceui->ui.powerRidge_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerRidge_l, AttType::powerRidge_l);
		connect(forceui->ui.powerArea_b,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerArea_b,AttType::powerArea_b);
		connect(forceui->ui.powerArea_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerArea_l,AttType::powerArea_l);
		connect(forceui->ui.powerArm_b,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerArm_b,AttType::powerArm_b);
		connect(forceui->ui.powerArm_l,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.powerArm_l,AttType::powerArm_l);

		connect(forceui->ui.biteRamp,SIGNAL(clicked()),attMap,SLOT(map()));		attMap->setMapping(forceui->ui.biteRamp, AttType::biteRamp);
		connect(forceui->ui.importObj,SIGNAL(clicked()),attMap,SLOT(map()));	attMap->setMapping(forceui->ui.importObj, AttType::importObj);

		connect(attMap, SIGNAL(mapped(int)), this, SLOT(slotBuildAtt(int)));
	}

	connect(forceui, SIGNAL(closing()),this->gla,SLOT(endEdit()));
	forceui->rebuildTree();
	forceui->show();

	if(!has_init)
	{
		CMN = teeth_tree[current_teeth]->nodeList[0];
		has_init = true;
	}
	
	curMode = ManuMode::None;
	inputValue = 0;
	original_Transform = CMN->m->cm.Tr;
	delta_Transform = vcg::Matrix44f::Identity();

	gla->update();
	rubberband.Reset();
	return true;
}

void ForceEditPlugin::EndEdit(MeshModel& m,GLArea* parent)
{
	if(forceui)
		delete forceui;
	forceui = 0;
	rubberband.Reset();
};


void ForceEditPlugin::slotChangeAttType()
{
	if (!CMN->isattachment)
	{
		return;
	}
	AttachmentNode *node = static_cast<AttachmentNode*>(CMN);
	int index = forceui->ui.cb_attLists->currentIndex();
	switch (index)
	{
	case 0:break;
	case 1:node->attType = AttType::rectVert;break;
	case 2:node->attType = AttType::rectHori;break;
	case 3:node->attType = AttType::rectInci;break;
	case 4:node->attType = AttType::rectGum;break;
	case 5:node->attType = AttType::rectMesio;break;
	case 6:node->attType = AttType::rectDistal;break;
	case 7:node->attType = AttType::ellipHori;break;
	case 8:node->attType = AttType::ellipVert;break;
	case 9:node->attType = AttType::optRotShort;break;
	case 10:node->attType = AttType::optRotLong;break;
	case 11:node->attType = AttType::optRootMesio;break;
	case 12:node->attType = AttType::optRootDistal;break;
	case 13:node->attType = AttType::optExt;break;
	case 14:node->attType = AttType::optDeep;break;
	case 15:node->attType = AttType::powerRidge_b;break;
	case 16:node->attType = AttType::powerRidge_l;break;
	case 17:node->attType = AttType::buttonCutout_l;break;
	case 18:node->attType = AttType::buttonCutout;break;
	case 19:node->attType = AttType::leftHook;break;
	case 20:node->attType = AttType::leftHook_l;break;
	case 21:node->attType = AttType::rightHook_l;break;
	case 22:node->attType = AttType::rightHook;break;
	case 23:node->attType = AttType::biteRamp;break;
	case 24:node->attType = AttType::buttonCutout_m;break;
	default:break;
	}
	return;
}


void trim_Brace(MeshModel* mm, const vector<Point3f>& vec_gumline, const int& which_jaw)
{
	// 遍历牙龈线上每个点，找出牙齿上离它最近的点 ********************************************************************  1
	vector<int> idx1;
	idx1.clear();

	for(int i = 0, gsize = vec_gumline.size(); i < gsize; i++)
	{
		int tmpi = -1;
		float mindis = 999.0f, dis = 0.0f;

		for(CMeshO::VertexIterator vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); ++vi)
		{
			if (!vi->IsD())
			{
				Point3f t = vi->P() - vec_gumline[i];
				if(fabs(t.X()) < 1.0f && fabs(t.Y()) < 1.0f && fabs(t.Z()) < 1.0f)	// 原来是3.0
				{
					dis = t.X()*t.X() + t.Y()*t.Y() + t.Z()*t.Z();	// 三个坐标平方和
					if(dis < mindis)
					{
						mindis = dis;
						tmpi = (int)(vi - mm->cm.vert.begin());
					}
				}
			}
		}
		if(tmpi == -1 || mindis > 0.2)	// 设置阈值，如果牙龈线的点离牙齿距离过远，则舍弃
			continue;
		//m.cm.vert[tmpi].C() = Color4b::Red;
		if(find(idx1.begin(), idx1.end(), tmpi) == idx1.end())
			idx1.push_back(tmpi);
	}
	int idx1len = idx1.size();
	//printf("idx1.size() = %d\n", idx1len);

	//for(int i = 0; i < idx1len; i++)
	//	mm->cm.vert[idx1[i]].C() = Color4b::Red;

	// 对idx1中每两个相邻点，用dijstra算法计算出最短路径，记录路径所有点的索引值 ***************************************** 3
	//printf("\n对idx1中每两个相邻点，用dijstra算法计算出最短路径，记录路径所有点的索引值\n");

	vector<int> idx2;
	idx2.clear();

	AbstractGraph* graph = new AbstractGraph(&(mm->cm));

	for(int i = 0; i < idx1len; i++)
	{
		vector<int> v;
		v.clear();
		int il = idx1[i];
		int ir = -1;
		if(i != idx1len - 1)
			ir = idx1[i + 1];
		else
			ir = idx1[0];

		vector<int>& neighbourList = graph->GetNeighbourList(il);
		if(find(neighbourList.begin(), neighbourList.end(), ir) != neighbourList.end())
		{
			//printf("两个点相连，直接 push_back 并继续\n");
			// il 肯定在 idx2中，所以只需要判断 ir
			if(find(idx2.begin(), idx2.end(), ir) == idx2.end())
				idx2.push_back(ir);
			continue;
		}

		GeodeticCalculator_Dijk* djk = new GeodeticCalculator_Dijk(*graph, il, ir);
		djk->Execute();
		v = djk->Get_Path();
		delete djk;


		for(int j = 0, vsize = v.size(); j < vsize; j++)
		{
			vector<int>::iterator ite = find(idx2.begin(), idx2.end(), v[j]);
			if(ite == idx2.end())
				idx2.push_back(v[j]);
		}
	}
	delete graph;

	//printf("idx2.size() = %d\n", idx2.size());

	tri::UpdateSelection<CMeshO>::ClearVertex(mm->cm);
	int idx2len = idx2.size();
	for(int i = 0; i < idx2len; i++)
	{
		//assert(!mm->cm.vert[idx2[i]].IsD());
		mm->cm.vert[idx2[i]].C() = Color4b::Red;
		mm->cm.vert[idx2[i]].SetS();
	}

	// 删除点和面	 *************************************************************************************************  4
	tri::UpdateSelection<CMeshO>::ClearFace(mm->cm);	// clear所有没有被delete的面
	//tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(m.cm);	// 有一个点被选中就将这个面setS()
	for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
	{
		if( !(*fi).IsD() && !(*fi).IsS())	
		{
			for(int i = 0; i < 3; i++)
			{
				if((*fi).V(i)->IsS())
				{
					fi->SetS();
					//fi->V(0)->C() = Color4b::Red;
					//fi->V(1)->C() = Color4b::Red;
					//fi->V(2)->C() = Color4b::Red;
					break;
				}
			}
		}
	}
	//for(int i = 0; i < idx1len; i++)
	//	mm->cm.vert[idx1[i]].C() = Color4b::Yellow;

	//for(int i = 0; i < idx2len; i++)
	//	mm->cm.vert[idx2[i]].C() = Color4b::Yellow;

	////// 从meshlab里拷过来的删除操作
	tri::UpdateSelection<CMeshO>::ClearVertex(mm->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(mm->cm);  
	for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
		if(!(*fi).IsD() && (*fi).IsS() )
			tri::Allocator<CMeshO>::DeleteFace(mm->cm, *fi);
	for(CMeshO::VertexIterator vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); ++vi)
		if(!(*vi).IsD() && (*vi).IsS() )
			tri::Allocator<CMeshO>::DeleteVertex(mm->cm, *vi);
	mm->clearDataMask(MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACEFLAGBORDER);
	////// 从meshlab里拷过来的删除操作


	// 计算出连通区域，利用种子面裁减掉牙龈部分 **********************************************************************  5
	tri::UpdateSelection<CMeshO>::ClearVertex(mm->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(mm->cm);

	mm->updateDataMask(MeshModel::MM_FACEFACETOPO);

	float maxZ = -999.0f;
	float minZ = 999.0f;
	int seed = -1;
	if(which_jaw == 0)
	{
		for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
		{
			if(!fi->IsD())
			{
				if(fi->V(1)->P().Z() > maxZ)
				{
					maxZ = fi->V(1)->P().Z();
					seed = (int)(fi - mm->cm.face.begin());
				}
			}
		}
	}
	else
	{
		for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
		{
			if(!fi->IsD())
			{
				if(fi->V(1)->P().Z() < minZ)
				{
					minZ = fi->V(1)->P().Z();
					seed = (int)(fi - mm->cm.face.begin());
				}
			}
		}
	}

	//assert(seed >= 0 && seed < mm->cm.face.size());
	CMeshO::FaceIterator ffi = mm->cm.face.begin() + seed;
	ffi->SetS();

	tri::UpdateSelection<CMeshO>::FaceConnectedFF(mm->cm);

	//////////////////////////////////////////// 删除连通区域

	///////////////////// 测试用，红色区域被扩充区域，即待裁剪区域
	//for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
	//	{
	//		if(!fi->IsD() && fi->IsS())
	//		{
	//			fi->V(0)->C() = Color4b::Red;
	//			fi->V(1)->C() = Color4b::Red;
	//			fi->V(2)->C() = Color4b::Red;
	//		}
	//	}
	//	mm->clearDataMask(MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACEFLAGBORDER);
	//	mm->updateDataMask(MeshModel::MM_FACEFACETOPO);
	/////////////////////// 测试用

	for(CMeshO::FaceIterator fi = mm->cm.face.begin(); fi != mm->cm.face.end(); ++fi)
	{
		if(!fi->IsD() && fi->IsS())
		{
			fi->V(0)->SetS();
			fi->V(1)->SetS();
			fi->V(2)->SetS();
			tri::Allocator<CMeshO>::DeleteFace(mm->cm, *fi);
		}
	}
	for(CMeshO::VertexIterator vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); ++vi)
	{
		if(!(*vi).IsD() && (*vi).IsS())
			tri::Allocator<CMeshO>::DeleteVertex(mm->cm, *vi);
	}
	mm->clearDataMask(MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACEFLAGBORDER);
	mm->updateDataMask(MeshModel::MM_FACEFACETOPO);

	////// BorderSmooth(mm);
	mm->updateDataMask(MeshModel::MM_ALL);//一定要放在前边

	tri::UpdateSelection<CMeshO>::FaceFromBorderFlag(mm->cm);
	tri::UpdateSelection<CMeshO>::VertexFromBorderFlag(mm->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mm->cm);
	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mm->cm);

	int stepSmoothNum = 10;
	float lambda = 1.0;
	float mu = -0.53;
	std::size_t cnt = tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(mm->cm);
	tri::Smooth<CMeshO>::VertexCoordTaubin(mm->cm,stepSmoothNum,lambda,mu,cnt>0);
	tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFace(mm->cm);	

	tri::UpdateBounding<CMeshO>::Box(mm->cm);
	tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mm->cm);
}

/////// 2018.5.25 起 - 计算所有期的带附件牙齿以及牙套，存入vector中 - GY
/////// 2018.5.25 起 - 计算所有期的带附件牙齿以及牙套，存入vector中 - GY
int ForceEditPlugin::compute_all_Teeth_and_Braces(int start, int end, int stride)
{
	long t1 = GetTickCount();

	// all_Teeth_and_Braces.clear(); 不释放内存，vecotr中如果存的是指针的话，不会调用对应析构函数，需要遍历，delete

	TeethTree* tt[2];
	tt[0] = teeth_tree[0];
	tt[1] = teeth_tree[1];

	int cnt = 0;
	for(int cur_step = start; cur_step <= end; cur_step += stride, ++cnt)//cur_step < end
	{
		printf("	———————————— 第 %d 步 开始 ————————————\n\n", cur_step);

		///////////////////////////////////////////// 1、计算下一期的牙齿、附件、牙龈
		printf("1、计算这一期的牙齿、附件、牙龈\n\n");

		for(int cti = 0; cti < 2; cti++)
		{
			if(cur_step >= tt[cti]->planstep)
				continue;
			tt[cti]->current_step = cur_step;

			// 更新牙齿以及附件
			foreach(ToothNode* nod, tt[cti]->seqtooth)
			{
				// 更新牙齿
				nod->nodetr = nod->move_state_matrix[cur_step];
				nod->m->cm.Tr = nod->nodetr;
				nod->UpdatefeaturePoint(nod->nodetr);
				nod->UpdateLocAxis(nod->nodetr);
				tri::UpdateBounding<CMeshO>::Box(nod->m->cm);
				tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(nod->m->cm);

				// 更新附件
				// 先默认一颗牙齿只有一个附件
				if(nod->attach.size() != 0 && nod->attach[0]->m->cm.vert.size() != 0)
				{
					nod->attach[0]->nodetr = nod->move_state_matrix[cur_step];
					nod->attach[0]->m->cm.Tr = nod->attach[0]->nodetr;
					nod->attach[0]->UpdatefeaturePoint(nod->attach[0]->nodetr);
					nod->attach[0]->UpdateLocAxis(nod->attach[0]->nodetr);
					tri::UpdateBounding<CMeshO>::Box(nod->attach[0]->m->cm);
					tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(nod->attach[0]->m->cm);
				}	
			}

			// 更新牙龈，主要为了得到里边的参数
			if(tt[cti]->hasgumfea)
			{
				foreach(ToothNode *mni,tt[cti]->nodeList)
				{
					if(mni->isgum)
						RTBGum::ExportGumCurvesandUnderCutSet(tt[cti],mni,tt[cti]->initseqtooth,cur_step);
				}
			}
		}

		////////////////////////////////// 2、读取gumcurvesorigin生成每一步的倒凹
		printf("2、读取生成每一步的倒凹\n\n");

		foreach(MeshModel* mm, md->meshList)
		{
			if(mm->label().contains("aligner"))
				md->delMesh(mm);
		}

		vector<Point3f> gumPoints;
		for(int cti = 0;cti < 2; cti++)
		{
			gumPoints.clear();

			foreach(ToothNode * nod, tt[cti]->nodeList)	// 一定是要nodeList，只有这里边有牙龈，seqtooth没有
			{
				if(!nod->isgum)
					continue;
				// std::vector<std::pair<std::string, std::vector<vcg::Point3f>> > undercutsets;// 填倒凹时候用的数据，每颗牙齿4个点

				if(cti == 0)
					printf("上颌有 %d 颗牙齿的牙龈点数据\n", nod->undercutsets.size());
				else
					printf("下颌有 %d 颗牙齿的牙龈点数据\n\n", nod->undercutsets.size());
				for(int i = 0; i < nod->undercutsets.size(); i++)
				{
					printf("牙龈点配置顺序!\n");
					//cout << "# "<< nod->undercutsets[i].first << endl;
					if(nod->undercutsets[i].second.size() != 0)
					{
						cout << nod->undercutsets[i].first << endl;
						vector<Point3f>& tmp = nod->undercutsets[i].second;
						assert(tmp.size() == 4);
						for(int j = 0; j < 4; j++)
							gumPoints.push_back(tmp[j]);
					}
				}
			}

			vector<MeshModel*> tooths;
			printf("牙齿顺序!\n");
			// 不知道为什么toothlab里和MQ里，下颌的顺序是反着的
			for(int i = 0; i < tt[cti]->seqtooth.size(); i++)
			{
				ToothNode* nod = tt[cti]->seqtooth[i];
				cout << nod->nodename << endl;
				MeshModel* m = nod->m;

				for(int j = 0; j < m->cm.vert.size(); j++)
					m->cm.vert[j].P() = nod->nodetr * m->cm.vert[j].P();

				tri::UpdateBounding<CMeshO>::Box(m->cm);
				tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(m->cm);
				tooths.push_back(m);
			}

			//if(cti == 0)	// U 的顺序是 UR7 到 UL7
			//{
			//	for(int i = 0; i < tt[cti]->seqtooth.size(); i++)
			//	{
			//		ToothNode* nod = tt[cti]->seqtooth[i];
			//		cout << nod->nodename << endl;
			//		MeshModel* m = nod->m;

			//		for(int j = 0; j < m->cm.vert.size(); j++)
			//			m->cm.vert[j].P() = nod->nodetr * m->cm.vert[j].P();

			//		tri::UpdateBounding<CMeshO>::Box(m->cm);
			//		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(m->cm);

			//		tooths.push_back(m);
			//	}
			//}
			//else	// L 的顺序是反着的，LL7 到 LR7
			//{
			//	//for(int i = tt[cti]->seqtooth.size() - 1; i >= 0; i--)
			//	{
			//		ToothNode* nod = tt[cti]->seqtooth[i];
			//		cout << nod->nodename << endl;
			//		MeshModel* m = nod->m;

			//		for(int j = 0; j < m->cm.vert.size(); j++)
			//			m->cm.vert[j].P() = nod->nodetr * m->cm.vert[j].P();

			//		tri::UpdateBounding<CMeshO>::Box(m->cm);
			//		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(m->cm);

			//		tooths.push_back(m);
			//	}
			//}

			vector<Point3f> readline;	// 没什么用，就是需要传一个？
			getUndetcut(tooths, gumPoints, md, readline, 0.4, 0.5);

			tooths.clear();

			///////////// 将点的坐标变回来
			if(cti == 0)
			{
				for(int i = 0; i < tt[cti]->seqtooth.size(); i++)
				{
					ToothNode* nod = tt[cti]->seqtooth[i];
					MeshModel* m = nod->m;

					vcg::Matrix44f re4;
					for(int i =0;i<3;i++)
						for(int j =0;j<3;j++)
							re4[i][j] = nod->nodetr[i][j];
					re4[0][3] = 0; re4[1][3] = 0; re4[2][3] = 0;
					re4[3][0] = 0; re4[3][1] = 0; re4[3][2] = 0; re4[3][3] = 1;

					nod->nodetr_rotation_transepose = re4;
					vcg::Transpose(nod->nodetr_rotation_transepose);

					Matrix44f inverse_translation;
					inverse_translation[0][0] = 1; inverse_translation[0][1] = 0; inverse_translation[0][2] = 0; inverse_translation[0][3] = -1 * nod->nodetr[0][3];
					inverse_translation[1][0] = 0; inverse_translation[1][1] = 1; inverse_translation[1][2] = 0; inverse_translation[1][3] = -1 * nod->nodetr[1][3];
					inverse_translation[2][0] = 0; inverse_translation[2][1] = 0; inverse_translation[2][2] = 1; inverse_translation[2][3] = -1 * nod->nodetr[2][3];
					inverse_translation[3][0] = 0; inverse_translation[3][1] = 0; inverse_translation[3][2] = 0; inverse_translation[3][3] = 1;

					for(int j = 0; j < m->cm.vert.size(); j++)
					{
						// 1、平移回坐标原点	// 2、旋转
						m->cm.vert[j].P() = nod->nodetr_rotation_transepose * (inverse_translation * m->cm.vert[j].P());
					}

					tri::UpdateBounding<CMeshO>::Box(m->cm);
					tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(m->cm);
				}
			}
			else
			{
				for(int i = tt[cti]->seqtooth.size() - 1; i >= 0; i--)
				{
					ToothNode* nod = tt[cti]->seqtooth[i];
					MeshModel* m = nod->m;

					vcg::Matrix44f re4;
					for(int i =0;i<3;i++)
						for(int j =0;j<3;j++)
							re4[i][j] = nod->nodetr[i][j];
					re4[0][3] = 0; re4[1][3] = 0; re4[2][3] = 0;
					re4[3][0] = 0; re4[3][1] = 0; re4[3][2] = 0; re4[3][3] = 1;

					nod->nodetr_rotation_transepose = re4;
					vcg::Transpose(nod->nodetr_rotation_transepose);

					Matrix44f inverse_translation;
					inverse_translation[0][0] = 1; inverse_translation[0][1] = 0; inverse_translation[0][2] = 0; inverse_translation[0][3] = -1 * nod->nodetr[0][3];
					inverse_translation[1][0] = 0; inverse_translation[1][1] = 1; inverse_translation[1][2] = 0; inverse_translation[1][3] = -1 * nod->nodetr[1][3];
					inverse_translation[2][0] = 0; inverse_translation[2][1] = 0; inverse_translation[2][2] = 1; inverse_translation[2][3] = -1 * nod->nodetr[2][3];
					inverse_translation[3][0] = 0; inverse_translation[3][1] = 0; inverse_translation[3][2] = 0; inverse_translation[3][3] = 1;

					for(int j = 0; j < m->cm.vert.size(); j++)
					{
						// 1、平移回坐标原点	// 2、旋转
						m->cm.vert[j].P() = nod->nodetr_rotation_transepose * (inverse_translation * m->cm.vert[j].P());
					}

					tri::UpdateBounding<CMeshO>::Box(m->cm);
					tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(m->cm);
				}
			}
			///////////// 将点的坐标变回来
		}

		///////////////////////////////////////////// 3、生成新的CSG数据，共4个，L_CSG_A、U_CSG_A为牙齿模型，L_CSG_B、U_CSG_B为倒凹以及附件模型
		printf("3、生成新的 CSG 数据\n");

		int s = 0;
		vector<Point3f> vertices;
		vector<Point3i> faces;
		vertices.clear();
		faces.clear();

		for(uint i = 0; i < tt[0]->seqtooth.size() - 1; i++)	// 到倒数第二个
		{
			cout << "计算" << tt[0]->seqtooth[i]->nodename << " 和 " << tt[0]->seqtooth[i + 1]->nodename << " 的碰撞" << endl;
			uint times = 0;
			while(twoTeethCollisionDetecte(tt[0]->seqtooth[i], tt[0]->seqtooth[i + 1]) == true && times < 15)	// 避免一直循环
			{
				++times;
				twoTeethCollisionRemove(tt[0]->seqtooth[i], tt[0]->seqtooth[i + 1]);
			}
			if(times == 15)
				printf("无法消除碰撞!!!!\n");
		}
		/////// 融合上颌所有牙齿[
		int cnt_UA = 0;
		foreach (ToothNode* nod, tt[0]->seqtooth)
		{
			++cnt_UA;
			s = vertices.size();
			CMeshO::VertexIterator vi;
			for(vi = nod->m->cm.vert.begin(); vi != nod->m->cm.vert.end(); ++vi)
				vertices.push_back(nod->nodetr * (vi->P()));

			CMeshO::FaceIterator fi;
			for (fi = nod->m->cm.face.begin(); fi != nod->m->cm.face.end();++fi)
			{
				CMeshO::FacePointer fp = &(*fi);
				Point3i a;
				a.X() = tri::Index(nod->m->cm, fp->V(0)) + s;
				a.Y() = tri::Index(nod->m->cm, fp->V(1)) + s;
				a.Z() = tri::Index(nod->m->cm, fp->V(2)) + s;
				faces.push_back(a);
			}
		}
		MeshModel* mUA = new MeshModel(md, "", "U_CSG_A.obj");
		vcg::tri::ConcaveCover<CMeshO>(mUA->cm, vertices, faces);
		printf("\t上颌牙齿融合完成，共 %d 颗牙齿，共 %d 个点，%d 个面!\n", cnt_UA, vertices.size(), faces.size());
		mUA->updateDataMask(MeshModel::MM_ALL);
		tri::UpdateBounding<CMeshO>::Box(mUA->cm);
		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mUA->cm);

		vertices.clear();
		faces.clear();
		s = 0;

		/////// 融合上颌所有倒凹以及附件
		int cnt_UB_1 = 0, cnt_UB_2 = 0;
		// 融合所有附件
		foreach (ToothNode* nod, tt[0]->seqtooth)
		{
			if(nod->attach.size() == 0 || nod->attach[0]->m->cm.vert.size() == 0)
				continue;

			++cnt_UB_1;
			s = vertices.size();
			CMeshO::VertexIterator vi;
			for(vi = nod->attach[0]->m->cm.vert.begin(); vi != nod->attach[0]->m->cm.vert.end(); ++vi)
				vertices.push_back(nod->attach[0]->nodetr * (vi->P()));

			CMeshO::FaceIterator fi;
			for (fi = nod->attach[0]->m->cm.face.begin(); fi != nod->attach[0]->m->cm.face.end();++fi)
			{
				CMeshO::FacePointer fp = &(*fi);
				Point3i a;
				a.X() = tri::Index(nod->attach[0]->m->cm, fp->V(0)) + s;
				a.Y() = tri::Index(nod->attach[0]->m->cm, fp->V(1)) + s;
				a.Z() = tri::Index(nod->attach[0]->m->cm, fp->V(2)) + s;
				faces.push_back(a);
			}
		}

		// 融合所有倒凹
		foreach(MeshModel* mm, md->meshList)
		{
			string label = string((const char*)mm->label().toLocal8Bit());
			if (label.find("-") != string::npos && (label.find("UL") != string::npos || label.find("UR") != string::npos))
			{
				++cnt_UB_2;

				s = vertices.size();
				CMeshO::VertexIterator vi;
				for(vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); ++vi)
					vertices.push_back(vi->P());

				CMeshO::FaceIterator fi;
				for (fi = mm->cm.face.begin(); fi != mm->cm.face.end();++fi)
				{
					CMeshO::FacePointer fp = &(*fi);
					Point3i a;
					a.X() = tri::Index(mm->cm,fp->V(0)) + s;
					a.Y() = tri::Index(mm->cm,fp->V(1)) + s;
					a.Z() = tri::Index(mm->cm,fp->V(2)) + s;
					faces.push_back(a);
				}
			}
		}
		MeshModel* mUB = new MeshModel(md, "", "U_CSG_B.obj");
		vcg::tri::ConcaveCover<CMeshO>(mUB->cm, vertices, faces);
		printf("\t上颌附件、倒凹融合完成，共 %d 个附件，%d 个倒凹，共 %d 个点，%d 个面!\n", cnt_UB_1, cnt_UB_2, vertices.size(), faces.size());
		mUB->updateDataMask(MeshModel::MM_ALL);
		tri::UpdateBounding<CMeshO>::Box(mUB->cm);
		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mUB->cm);

		vertices.clear();
		faces.clear();
		s = 0;

		/////// 融合下颌所有牙齿
		int cnt_LA = 0;
		for(uint i = 0; i < tt[1]->seqtooth.size() - 1; i++)	// 到倒数第二个
		{
			cout << "计算" << tt[1]->seqtooth[i]->nodename << " 和 " << tt[1]->seqtooth[i + 1]->nodename << " 的碰撞" << endl;
			uint times = 0;
			while(twoTeethCollisionDetecte(tt[1]->seqtooth[i], tt[1]->seqtooth[i + 1]) == true && times < 15)// && times < 8
			{
				++times;
				twoTeethCollisionRemove(tt[1]->seqtooth[i], tt[1]->seqtooth[i + 1]);
			}
			if(times == 15)
				printf("无法消除碰撞!!!!\n");
		}

		foreach (ToothNode* nod, tt[1]->seqtooth)
		{
			++cnt_LA;
			s = vertices.size();
			CMeshO::VertexIterator vi;
			for(vi = nod->m->cm.vert.begin(); vi != nod->m->cm.vert.end(); ++vi)
				vertices.push_back(nod->nodetr * (vi->P()));

			CMeshO::FaceIterator fi;
			for (fi = nod->m->cm.face.begin(); fi != nod->m->cm.face.end();++fi)
			{
				CMeshO::FacePointer fp = &(*fi);
				Point3i a;
				a.X() = tri::Index(nod->m->cm, fp->V(0)) + s;
				a.Y() = tri::Index(nod->m->cm, fp->V(1)) + s;
				a.Z() = tri::Index(nod->m->cm, fp->V(2)) + s;
				faces.push_back(a);
			}
		}
		MeshModel* mLA = new MeshModel(md, "", "L_CSG_A.obj");
		vcg::tri::ConcaveCover<CMeshO>(mLA->cm, vertices, faces);
		printf("\t下颌牙齿融合完成，共 %d 颗牙齿，共 %d 个点，%d 个面!\n", cnt_LA, vertices.size(), faces.size());
		mUA->updateDataMask(MeshModel::MM_ALL);
		tri::UpdateBounding<CMeshO>::Box(mLA->cm);
		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mLA->cm);

		vertices.clear();
		faces.clear();
		s = 0;

		/////// 融合下颌所有倒凹以及附件
		int cnt_LB_1 = 0, cnt_LB_2 = 0;
		// 融合所有附件
		foreach (ToothNode* nod, tt[1]->seqtooth)
		{
			if(nod->attach.size() == 0 || nod->attach[0]->m->cm.vert.size() == 0)
				continue;

			++cnt_LB_1;
			s = vertices.size();
			CMeshO::VertexIterator vi;
			for(vi = nod->attach[0]->m->cm.vert.begin(); vi != nod->attach[0]->m->cm.vert.end(); ++vi)
				vertices.push_back(nod->attach[0]->nodetr * (vi->P()));

			CMeshO::FaceIterator fi;
			for (fi = nod->attach[0]->m->cm.face.begin(); fi != nod->attach[0]->m->cm.face.end();++fi)
			{
				CMeshO::FacePointer fp = &(*fi);
				Point3i a;
				a.X() = tri::Index(nod->attach[0]->m->cm, fp->V(0)) + s;
				a.Y() = tri::Index(nod->attach[0]->m->cm, fp->V(1)) + s;
				a.Z() = tri::Index(nod->attach[0]->m->cm, fp->V(2)) + s;
				faces.push_back(a);
			}
		}

		// 融合所有倒凹
		foreach(MeshModel* mm, md->meshList)
		{
			string label = string((const char*)mm->label().toLocal8Bit());
			if (label.find("-") != string::npos && (label.find("LL") != string::npos || label.find("LR") != string::npos))
			{
				++cnt_LB_2;

				s = vertices.size();
				CMeshO::VertexIterator vi;
				for(vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); ++vi)
					vertices.push_back(vi->P());

				CMeshO::FaceIterator fi;
				for (fi = mm->cm.face.begin(); fi != mm->cm.face.end();++fi)
				{
					CMeshO::FacePointer fp = &(*fi);
					Point3i a;
					a.X() = tri::Index(mm->cm,fp->V(0)) + s;
					a.Y() = tri::Index(mm->cm,fp->V(1)) + s;
					a.Z() = tri::Index(mm->cm,fp->V(2)) + s;
					faces.push_back(a);
				}
			}
		}
		MeshModel* mLB = new MeshModel(md, "", "L_CSG_B.obj");
		vcg::tri::ConcaveCover<CMeshO>(mLB->cm, vertices, faces);
		printf("\t上颌附件、倒凹融合完成，共 %d 个附件，%d 个倒凹，共 %d 个点，%d 个面!\n", cnt_LB_1, cnt_LB_2, vertices.size(), faces.size());
		mLB->updateDataMask(MeshModel::MM_ALL);
		tri::UpdateBounding<CMeshO>::Box(mLB->cm);
		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mLB->cm);

		vertices.clear();
		faces.clear();

		printf("成功生成四个 CSG 文件！\n\n");

		///////////////////////////////////////////// 4、对上下颌CSG模型进行 CSG Union 操作
		printf("4、对上下颌 CSG 模型进行 CSG Union 操作\n");

		CSG_Union(mUA, mUB, md, 0);
		CSG_Union(mLA, mLB, md, 1);

		// CSG 结束了，生成了两个牙套，4个CSG文件以及倒凹就都没有用了
		delete mUA;
		delete mUB;
		delete mLA;
		delete mLB;

		// 将倒凹删除，这样剩下的带有LR等的文件都是牙齿
		foreach(MeshModel* mm, md->meshList)
		{
			if(mm->label().contains("-"))
				md->delMesh(mm);
		}

		printf("\tCSG结束!\n\n");

		///////////////////////////////////////////// 5、对CSG 结果（两个aligner模型），进行裁剪（根据牙龈的 gumcurvesorign 属性）
		printf("5、对 CSG 结果进行裁剪\n");
		// 在程序中得到牙龈点数据

		vector<Point3f> gumpoints_L, gumpoints_U;

		gumpoints_L.clear();
		gumpoints_U.clear();
		foreach(ToothNode* nod, tt[0]->nodeList)
		{
			if(nod->isgum && nod->gumcurvesorign.size() != 0)
			{
				gumpoints_U = nod->gumcurvesorign;
				break;
			}
		}
		foreach(ToothNode* nod, tt[1]->nodeList)
		{
			if(nod->isgum && nod->gumcurvesorign.size() != 0)
			{
				gumpoints_L = nod->gumcurvesorign;
				break;
			}
		}
		assert(gumpoints_L.size() != 0 && gumpoints_U.size() != 0);
		//printf("下颌有 %d 个牙龈点， 上颌有 %d 个牙龈点\n", gumpoints_L.size(), gumpoints_U.size());

		int cnt1 = 0;
		foreach(MeshModel* mm, md->meshList)
		{
			if(mm->label().contains("U_aligner"))
			{
				mm->cm.vn = mm->cm.vert.size();
				mm->cm.fn = mm->cm.face.size();
				printf("\t上颌牙套裁剪，原有 %d 个点，%d 个面\n", mm->cm.vn, mm->cm.fn);

				trim_Brace(mm, gumpoints_U, 0);
				printf("\t上颌裁剪完成，剩余 %d 个点，%d 个面\n", mm->cm.vn, mm->cm.fn);
				++cnt1;
			}
			else if(mm->label().contains("L_aligner"))
			{
				mm->cm.vn = mm->cm.vert.size();
				mm->cm.fn = mm->cm.face.size();
				printf("\t下颌牙套裁剪，原有 %d 个点，%d 个面\n", mm->cm.vn, mm->cm.fn);

				trim_Brace(mm, gumpoints_L, 1);
				printf("\t下颌裁剪完成，剩余 %d 个点，%d 个面\n", mm->cm.vn, mm->cm.fn);
				++cnt1;
			}
			if(cnt1 == 2)
				break;
		}

		printf("	牙套裁剪完成!\n\n");

		///////////////////////////////////////////// 6、将牙套存入vector中
		printf("6、将牙套存入vector中\n");

		int mm_id = 0;

		ToothNode* aligner_U;
		ToothNode* aligner_L;
		foreach(MeshModel* mm, md->meshList)
		{
			if(mm->label().contains("aligner"))
			{
				/////////////////// 创建一个新的 meshmodel实体，用它创建 ToothNode 指针

				// MeshModel * MeshDocument::addNewMesh(QString fullPath, QString label, bool setAsCurrent)
				MeshModel* mm_tmp = new MeshModel(md, mm->fullName(), mm->label());

				// vcglib\vcg\complex\trimesh 这个文件夹下的 append.h
				// static void Mesh(MeshLeft& ml, MeshRight& mr, const bool selected = false)
				// 对自己创建的模型不好使，不知道什么原因，在copy faces时候出错，导入的模型都是好使的，所以自己写了一个 Copy_CMeshO 方法
				vcg::tri::Append<CMeshO, CMeshO>::Copy_CMeshO(mm_tmp->cm, mm->cm);

				/////////////////// 创建一个新的 meshmodel实体，用它创建 ToothNode 指针

				// nod->nodename 形式 LR5, mm->label().toStdString() 形式 LR5.obj
				string tmp_nodename = string((const char*)mm_tmp->label().toLocal8Bit());
				tmp_nodename = tmp_nodename.substr(0, tmp_nodename.length() - 4);

				// ToothNode(MeshModel *_m, int _id,std::string _nodename)
				if(mm->label().contains("U_"))
				{
					aligner_U = new ToothNode(mm_tmp, mm_id++, tmp_nodename);
					aligner_U->isalgner = true;
				}
				else if(mm->label().contains("L"))
				{
					aligner_L = new ToothNode(mm_tmp, mm_id++, tmp_nodename);
					aligner_L->isalgner = true;
				}
			}
		}
		vector<ToothNode* > current_step_aligners;

		current_step_aligners.push_back(aligner_U);
		printf("\tU_aligner 存入vector，共 %d 个点，%d 个面\n", aligner_U->m->cm.vert.size(), aligner_U->m->cm.face.size());
		current_step_aligners.push_back(aligner_L);
		printf("\tL_aligner 存入vector，共 %d 个点，%d 个面\n", aligner_L->m->cm.vert.size(), aligner_L->m->cm.face.size());

		all_Braces.push_back(current_step_aligners);

		///////////////////////////////////////////// 这一步完成
		printf("\n	———————————— 第 %d 步 结束 ————————————\n\n", cur_step);
		

	}

	///////////////////////////////////////////// 将牙齿和附件的变换矩阵变回最开始的
	for(int cti = 0; cti < 2;cti++)
	{	
		tt[cti]->current_step = 0;

		// 还原牙齿以及附件（回到最开始，第0期）
		foreach(ToothNode* nod, tt[cti]->seqtooth)
		{
			// 还原牙齿
			nod->nodetr = nod->move_state_matrix[0];
			nod->m->cm.Tr = nod->nodetr;
			nod->UpdatefeaturePoint(nod->nodetr);
			nod->UpdateLocAxis(nod->nodetr);
			tri::UpdateBounding<CMeshO>::Box(nod->m->cm);
			tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(nod->m->cm);

			// 还原附件
			if(nod->attach.size() != 0)
			{
				nod->attach[0]->nodetr = nod->move_state_matrix[0];
				nod->attach[0]->m->cm.Tr = nod->attach[0]->nodetr;
				nod->attach[0]->UpdatefeaturePoint(nod->attach[0]->nodetr);
				nod->attach[0]->UpdateLocAxis(nod->attach[0]->nodetr);
				tri::UpdateBounding<CMeshO>::Box(nod->attach[0]->m->cm);
				tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(nod->attach[0]->m->cm);
			}	
		}

		// 没有还原牙龈 2018.6.5
	}

	for(int i = 0; i < all_Braces.size(); i++)
	{
		printf("第 %d 期的 all_Braces\n", i + 1);	// 最开始就是第一期的，第0期的牙套没意义

		MeshModel* tmm = all_Braces[i][0]->m;
		printf("\t上颌牙套有 %d 个点，%d 个面\n", tmm->cm.vert.size(), tmm->cm.face.size());

		tmm = all_Braces[i][1]->m;
		printf("\t下颌牙套有 %d 个点，%d 个面\n", tmm->cm.vert.size(), tmm->cm.face.size());

		printf("\n");
	}

	long t2 = GetTickCount();
	long time_used = (t2 - t1) * 1.0 / 1000;
	int minutes = time_used / 60;
	int seconds = time_used % 60;
	printf("\n上下颌牙套计算完成，用时: %dm %ds\n\n", minutes, seconds);

	return cnt;
}

bool ForceEditPlugin::twoTeethCollisionDetecte(ToothNode* tn1, ToothNode* tn2)
{
	bool collision1 = false;
	bool collision2 = false;

	std::vector<Collide_Result* > res1, res2;

	ConlisionSimulate* cp1 = new ConlisionSimulate(tn1, tn2, tn1->nodetr, tn2->nodetr);//mv1 get red
	collision1 = cp1->DetectCollision(res1);
	delete cp1;

	ConlisionSimulate* cp2 = new ConlisionSimulate(tn2, tn1, tn2->nodetr, tn1->nodetr);//mv1 get red	
	collision2 = cp2->DetectCollision(res2);
	delete cp2;

	if (collision1 || collision2)
	{
		cout << "\t" << tn1->nodename << " 和 " << tn2->nodename << " 发生碰撞" << endl;
		return true;
	}
	else
		return false;
}

void ForceEditPlugin::twoTeethCollisionRemove(ToothNode* tn1, ToothNode* tn2)
{
	MeshModel *mm_smooth;
	CMeshO::VertexIterator vi;

	mm_smooth = tn1->m;
	tri::UpdateSelection<CMeshO>::ClearVertex(mm_smooth->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(mm_smooth->cm);

	vi = mm_smooth->cm.vert.begin();
	for (int i = 0; i <mm_smooth->cm.vert.size();i++)
	{
		vcg::Color4f c = Color4f::Construct(vi[i].C());
		if (c == Color4<float>(Color4<float>::Red))
		{
			vi[i].SetS();
		}
	}

	tri::UpdateFlags<CMeshO>::FaceBorderFromNone(mm_smooth->cm);

	// 扩充一圈set的区域
	for(CMeshO::FaceIterator fi = mm_smooth->cm.face.begin(); fi != mm_smooth->cm.face.end(); fi++)
	{
		if((*fi).V(0)->IsS() || (*fi).V(1)->IsS() || (*fi).V(2)->IsS())
			fi->SetS();
	}
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mm_smooth->cm);
	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mm_smooth->cm);

	tri::Smooth<CMeshO>::VertexCoordLaplacian(mm_smooth->cm, 5, true);
	tri::UpdateBounding<CMeshO>::Box(mm_smooth->cm);
	if(mm_smooth->cm.fn > 0) 
	{
		tri::UpdateNormals<CMeshO>::PerFaceNormalized(mm_smooth->cm);
		tri::UpdateNormals<CMeshO>::PerVertexAngleWeighted(mm_smooth->cm);
	}

	mm_smooth = tn2->m;
	tri::UpdateSelection<CMeshO>::ClearVertex(mm_smooth->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(mm_smooth->cm);

	vi = mm_smooth->cm.vert.begin();
	for (int i = 0; i <mm_smooth->cm.vert.size();i++)
	{
		vcg::Color4f c = Color4f::Construct(vi[i].C());
		if (c == Color4<float>(Color4<float>::Red))
		{
			vi[i].SetS();
		}
	}

	tri::UpdateFlags<CMeshO>::FaceBorderFromNone(mm_smooth->cm);

	// 扩充一圈set的区域
	for(CMeshO::FaceIterator fi = mm_smooth->cm.face.begin(); fi != mm_smooth->cm.face.end(); fi++)
	{
		if((*fi).V(0)->IsS() || (*fi).V(1)->IsS() || (*fi).V(2)->IsS())
			fi->SetS();
	}
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mm_smooth->cm);
	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mm_smooth->cm);

	tri::Smooth<CMeshO>::VertexCoordLaplacian(mm_smooth->cm, 5, true);
	tri::UpdateBounding<CMeshO>::Box(mm_smooth->cm);
	if(mm_smooth->cm.fn > 0) 
	{
		tri::UpdateNormals<CMeshO>::PerFaceNormalized(mm_smooth->cm);
		tri::UpdateNormals<CMeshO>::PerVertexAngleWeighted(mm_smooth->cm);
	}
}


void ForceEditPlugin::slotAttachOptim()	// A5E附件优化，现在变成不进行优化，只输出一次的误差，2018.8.20
{
	int start = 1;
	int end = std::max(teeth_tree[0]->planstep, teeth_tree[1]->planstep) / 2;
	int stride = 1;


	printf("开始计算误差!\n");
	long t1 = GetTickCount();

	map<string, vector<float>> error_map;
	//typedef pair<string, vector<int>> error_typr;
	// 上颌牙齿
	error_map.insert(pair<string, vector<float>>("UR7", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR6", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR5", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR4", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR3", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR2", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UR1", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL7", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL6", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL5", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL4", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL3", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL2", vector<float>()));
	error_map.insert(pair<string, vector<float>>("UL1", vector<float>()));

	// 下颌牙齿
	error_map.insert(pair<string, vector<float>>("LR7", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR6", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR5", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR4", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR3", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR2", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LR1", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL7", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL6", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL5", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL4", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL3", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL2", vector<float>()));
	error_map.insert(pair<string, vector<float>>("LL1", vector<float>()));

	// 计算每一期上下颌牙套，存入all_Braces里

	//for(int i = 0; i < 2; i++)	// 优化次数，优化几次应该有几个新附件
	//{
	// 每次附件优化（大循环）需要将 all_Braces 清空
	//for(int m = 0; m < all_Braces.size(); m++)
	//{
	//	if(all_Braces[m].size() != 2)
	//		printf("\n\n牙套数量出错!\n");
	//	for(int n = 0; n < 2; n++)
	//		delete all_Braces[m][n];
	//}
	//all_Braces.clear();

	printf("\n计算每一期上下颌牙套，存入all_Braces里!\n\n");
	int steps = compute_all_Teeth_and_Braces(start, end, stride);	// GY

	printf("牙套计算完成，开始进行附件优化!\n\n");
	compute_all_Force_and_Torque(error_map, steps);

	//	printf("\n\n第 %d 次附件优化完成!!!\n\n", i);
	//}

	ofstream error_file;
	//  error_file.open("../error_file.txt", ios::out | ios::app);	// 从文件末尾开始写
	error_file.open("../error_file.txt");	// 覆盖之前的

	map< string, vector<float> >::iterator it;  
	for(it = error_map.begin(); it != error_map.end(); it++)
	{
		error_file << endl << it->first << endl;	// 牙齿名称
		error_file << "error1 - 力的误差";
		for(int i = 0; i < it->second.size(); i = i + 2)
			error_file << " " << it->second[i];
		error_file << endl << "error2 - 力矩的误差";
		for(int i = 1; i < it->second.size(); i = i + 2)
			error_file << " " << it->second[i];
		error_file << endl;
	}

	error_file.close();

	long t2 = GetTickCount();
	long time_used = (t2 - t1) * 1.0 / 1000;
	int minutes = time_used / 60;
	int seconds = time_used % 60;
	printf("**************************************\n");
	printf("\n误差计算完成，输出为error_file.txt文件，用时: %dm %ds\n\n", minutes, seconds);
	//printf("\n共进行 %d 次附件优化，用时: %dm %ds\n\n", i, minutes, seconds);
	printf("**************************************\n");
	
}

int itn = 0;
void ForceEditPlugin::compute_all_Force_and_Torque(map<string, vector<float>>& error_map, int steps)
{
	for(int cur_teeth = 0; cur_teeth <= 1; cur_teeth++)
	{
		for (size_t i = 0; i < steps; i++)// max_step - 1
		{
			printf("计算第 %d 期的牙齿受第 %d 期牙套的力和力矩, current_teeth = %d\n", i, i + 1, cur_teeth);
			slotComputeForceofsteps(i, cur_teeth);//在中间修改当前teethtree中的牙套对其所受的力和力矩
		}
		//误差判别
		itn++;

		foreach(ToothNode *tni,teeth_tree[cur_teeth]->seqtooth)
		{
			//printf("\n%s\n", (char*)(tni->nodename.c_str()));
			Point3f n_TF_Normal = tni->n_TargetForce.Normalize();
			Point3f AAF_Normal = tni->AllAlignerForce.Normalize();
			Point3f n_TT_Normal = tni->n_TargetTorque.Normalize();
			Point3f AAT_Normal = tni->AllAlignerTorque.Normalize();

			//printf("\tn_TargetForce = (%f, %f, %f)\n", n_TF_Normal.X(), n_TF_Normal.Y(), n_TF_Normal.Z());
			//printf("\tAllAlignerForce = (%f, %f, %f)\n", AAF_Normal.X(), AAF_Normal.Y(), AAF_Normal.Z());
			//printf("\tn_TargetTorque = (%f, %f, %f)\n", n_TT_Normal.X(), n_TT_Normal.Y(), n_TT_Normal.Z());
			//printf("\tAllAlignerTorque = (%f, %f, %f)\n\n", AAT_Normal.X(), AAT_Normal.Y(), AAT_Normal.Z());

			//printf("\tn_TF_Normal * AAF_Normal = %f\n", n_TF_Normal * AAF_Normal);
			//printf("\tacos(n_TF_Normal * AAF_Normal) = %lf\n", acos(n_TF_Normal * AAF_Normal));

			float error1 = (tni->n_TargetForce).Norm() == 0 ? 0 : acos(n_TF_Normal * AAF_Normal);
			float error2 = (tni->n_TargetTorque).Norm() == 0 ? 0 : acos(n_TT_Normal * AAT_Normal);

			error_map[tni->nodename].push_back(error1);
			error_map[tni->nodename].push_back(error2);

			if (abs(error1) + abs(error2) > 3.1415 / 3.0)
			{
				//重新计算该ToothNode的附件
				cout<<tni->nodename<<endl;
				ComputeAttachForOneTooth(tni);
				printf("所有力的误差 = %f\n", error1);
				printf("所有力矩的误差 = %f\n\n", error2);
				//printf("=========================\n");
			}
		}
	}
}


// 生成第 i 期带附件牙齿，与第 i+1 期牙套进行受力
void ForceEditPlugin::slotComputeForceofsteps(int i, int cur_teeth)
{
	//计算牙齿的力

	std::vector<ToothNode *> aligner;
	aligner.clear();
	aligner.push_back(all_Braces[i][cur_teeth]);//读入aligner文件

	// 生成第 i 期带附件牙齿 - 开始
	vector<ToothNode* > cur_tooth;
	foreach(ToothNode* nod, teeth_tree[cur_teeth]->seqtooth)
	{
		MeshModel* tmp_m = new MeshModel(md, "", nod->m->label());
		// 如果想要创建新的ToothNode，必须先创建MeshModel
		// ToothNode(MeshModel *_m, int _id,std::string _nodename)
		ToothNode* tmp = new ToothNode(tmp_m, 0, nod->nodename);

		// 如果有对应附件且附件不是空的，则CSG融合，并赋予tmp
		if(nod->attach.size() != 0 && nod->attach[0]->m->cm.vert.size() != 0)
		{
			// CSG融合牙齿以及其对应附件，目前认为一颗牙齿最多有一个附件
			//printf("\t融合前 %s 有 %d 个点，附件有 %d 个点\n", (char*)(nod->nodename.c_str()), nod->m->cm.vert.size(), nod->attach[0]->m->cm.vert.size());

			CSG_Union_Teeth_and_Att_matrix(nod->m, nod->attach[0]->m, nod->move_state_matrix[cur_teeth], tmp->m, cur_teeth);
		}
		else
		{
			vcg::tri::Append<CMeshO, CMeshO>::Mesh(tmp->m->cm, nod->m->cm);
		}

		tmp->m->cm.vn = tmp->m->cm.vert.size();
		tmp->m->cm.fn = tmp->m->cm.face.size();
		// printf("\t有 %d 个点，%d 个面\n\n", tmp->m->cm.vn, tmp->m->cm.fn);

		tmp->m->updateDataMask(MeshModel::MM_ALL);
		tri::UpdateBounding<CMeshO>::Box(tmp->m->cm);
		tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(tmp->m->cm);

		cur_tooth.push_back(tmp);
	}
	printf("第 %d 期 %d 颗牙齿存入cur_tooth中用于碰撞检测\n\n", i, cur_tooth.size());

	// 生成第 i 期带附件牙齿 - 结束

	std::map<CVertexO *,double> res;
	//CMN->m->updateDataMask(MeshModel::MM_VERTMARK);

	std::vector<ToothNode *> teeths;	// 其实每次只放一个，下次使用时先clear
	foreach(ToothNode* tni, cur_tooth)
	{
		teeths.clear();

		//cout << "计算 " << tni->nodename << "与牙套的碰撞" << endl;
		teeths.push_back(tni);//读入牙齿和牙龈
		ConlisionSimulate *csi = new ConlisionSimulate(aligner, teeths);//将牙套与牙齿进行碰撞检测

		SetCurrentNode(tni);

		CMN->m->updateDataMask(MeshModel::MM_VERTMARK);
		CMN->m->updateDataMask(MeshModel::MM_VERTFACETOPO);
		CMN->m->updateDataMask(MeshModel::MM_FACEFACETOPO);

		res.clear();

		csi->Collidevalue(res);//将res进行赋值

		ComputeAlignerForceofTooth(res);

		//需要在原来的牙齿中，将力和力矩记录下来

		if(teeth_tree[cur_teeth]->GetToothByName(CMN->nodename) == false)
			printf("Wrong Name!\n");
		//printf("AllAlignerForce = (%f, %f, %f)\n", tmp.X(), tmp.Y(), tmp.Z());
		teeth_tree[cur_teeth]->GetToothByName(CMN->nodename)->AllAlignerForce += CMN->AlignerForce;
		teeth_tree[cur_teeth]->GetToothByName(CMN->nodename)->AllAlignerTorque += CMN->AlignerTorque;

		delete csi;
	}
	// 现在ToothNode的析构函数还是没有，所以先不释放空间 2018.6.13 - GY
	foreach(ToothNode* tni, cur_tooth)
	{
		delete tni->m;
		tni = NULL;
	}
	cur_tooth.clear();
}

void ForceEditPlugin::ComputeAlignerForceofTooth(std::map<CVertexO *,double> &res)
{
	CMeshO::VertexIterator vi;
	for(vi = CMN->m->cm.vert.begin(); vi != CMN->m->cm.vert.end(); ++vi)
		if( !(*vi).IsD() )
			Mark(&(*vi),vcg::ForceMark::FMN);
	std::map<CVertexO *,double>::const_iterator mci;
	vcg::Point3f temp;
	for(mci = res.begin();mci!= res.end();mci++)
	{
		CVertexO *cvi = (*mci).first;
		Mark(cvi,vcg::ForceMark::FMY);

		temp=cvi->P();
		CMN->dot.push_back(temp);
	}
	CMeshO::FaceIterator fi;
	std::vector<ForceSet> force_set;
	force_set.clear();
	for(fi = CMN->m->cm.face.begin(); fi != CMN->m->cm.face.end(); ++fi)
	{
		int num =0;
		ForceSet forcei;//这个数据结构有点迷
		forcei.f_beginpos.SetZero();
		forcei.f_direction.SetZero();
		forcei.f_value = 0.0;
		double facesurface = ((((*fi).V(1)->P() - (*fi).V(0)->P()) ^ ((*fi).V(2)->P() - (*fi).V(0)->P())).Norm())/2;//面片的面积
		double avd=0;
		for(int i =0; i< 3;i++)
		{				
			CVertexO *cvi = (*fi).V(i);
			if(cvi->IMark() == vcg::ForceMark::FMY)
			{				
				//if(res.count(cvi)){
				avd+=res[cvi];
				//cout<<"碰撞深度："<<res[cvi]<<endl;
				num++;
				//cout<<"Holyshit!!!!!"<<endl;
			}
			forcei.f_beginpos += cvi->P();
		}
		if(num >=2)//同一个面片被标记2次才算碰撞到
		{
			avd/=num;//碰撞的平均深度
			forcei.f_direction = -((*fi).N());//力的方向？
			forcei.f_value = Kc*avd*facesurface;
			//forcei.f_direction.Normalize();去除标准化步骤，向量的模直接表示大小
			forcei.f_beginpos /= 3;//面片的重心
			force_set.push_back(forcei);
		}
	}
	//开始计算每一个面片的力和力矩
	ForceVector.SetZero();//记录所有的碰撞点的力
	Point3f Torquevector;
	Torquevector.SetZero();//记录所有碰撞点的力矩
	for(int i =0;i< force_set.size();i++)
	{
		ForceSet forcei = force_set[i];
		Point3f f = forcei.f_direction*forcei.f_value;//面片的力
		ForceVector += f;
		Torquevector += (forcei.f_beginpos-CMN->barycentric_coord) ^ f;//面片的力矩
	}
	CMN->AlignerForce = ForceVector;
	CMN->AlignerTorque = -Torquevector;
	//cout << "ForceVector.Norm() = " << ForceVector.Norm() <<endl;
	//cout << "Torquevector.Norm() = " << Torquevector.Norm() << endl;
	return;
}


void ForceEditPlugin::ComputeAttachForOneTooth(ToothNode* tni){
	tni->attach.clear();
	for (int j = 0; j < tni->attloc.size(); j++)
	{
		int count = 0;
		while (teeth_tree[current_teeth]->GetAttByName(tni->qnodename + QString::number(count, 10)))
		{
			count++;
		}
		QString newNodeName = tni->qnodename + QString::number(count, 10);
		MeshModel* m = md->addNewMesh(newNodeName, false);
		AttachmentNode *newNode = new AttachmentNode(m, teeth_tree[current_teeth]->attNodes.size(), newNodeName, tni);
		std::string s = tni->nodename;
		Matrix44f ro; ro.SetIdentity();
		newNode->setPickedInfo(tni->attloc[j], tni->attforces[j]);
		Point3f vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1), Ou(0, 0, 0), p = tni->attloc[j], vn(0, 0, 1), vvx(1, 0, 0), vl(0, 0, 0);
		float pra1 = 0, pra2 = 0, alpha = 0.9;
		if (s == "UL6" || s == "UL7" || s == "UR6" || s == "UR7")
		{
			if (tni->toothFeature.crownHeight > 3 && teeth_tree[current_teeth]->hasanchorage)
			{
				if (abs(tni->y_tran) > 2) 
				{
					newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->extint) > 2 || abs(tni->rotation_theta) > 20)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->x_tran) > 2 || abs(tni->x_rotataion) > 20 || abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					vvx = vn^vvx;
				}
				else
				{
					newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot, tni->pca[1]);;
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->fb;
					pra1 = 0; pra2 = 0;
				}
			}
		}
		else if (s == "LL6" || s == "LL7" || s == "LR6" || s == "LR7")
		{
			if (tni->toothFeature.crownHeight > 3 && teeth_tree[current_teeth]->hasanchorage)
			{
				if (abs(tni->y_tran) > 2) 
				{
					newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->extint) > 2 || abs(tni->x_tran) > 2 || abs(tni->y_rotataion) > 20)
				{
					newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
				else if (abs(tni->x_rotataion) > 20 || abs(tni->rotation_theta) > 20)
				{
					newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
					pra1 = 70; pra2 = 0;
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					vvx = vn^vvx;
				}
				else {
					newNode->setType(AttType::ellipHori);//垂直矩形附件,根据附着位置的法向量进行粘贴
					ro = ArbRot(FPi / 2.0, Ou, vy);
					vn = NearestVer(tni, tni->toothFeature.lateralPivot);
					Point3f vv = vx.Normalize() ^ vn.Normalize();
					float rotvalue = acos(vx.Normalize()*vn.Normalize());
					ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
					vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
					p = tni->toothFeature.lateralPivot;
					pra1 = 0; pra2 = 0;
				}
			}
		}
		else if (s == "UL4" || s == "UL5" || s == "UR4" || s == "UR5" || s == "LL4" || s == "LL5" || s == "LR4" || s == "LR5")
		{
			if (abs(tni->extint) > 2 || abs(tni->y_tran) > 2)
			{
				newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (abs(tni->x_tran) > 2)
			{
				newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
				pra1 = 70; pra2 = 0;
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				vvx = vn^vvx;
			}
			else if (abs(tni->x_rotataion) > 20)
			{
				newNode->setType(AttType::optRoot);//优化控根
				vn = NearestVer(tni, tni->attloc[j]);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv, 0.7);//绕着经过p的vv旋转rotvalue(弧度)

				p = tni->attloc[j];
				pra1 = 0; pra2 = 0;
				if ((p - tni->barycentric_coord)*(tni->pca[2]) < 0)
				{
					vvx = tni->pca[2] - vn / (1.0 / (tni->pca[2] * vn));
				}
				else {
					vvx = -(tni->pca[2] - vn / (1.0 / (tni->pca[2] * vn)));
				}

			}
			else if (abs(tni->rotation_theta) > 20) {
				newNode->setType(AttType::optRotLong);//优化旋转
				vn = NearestVer(tni, tni->attloc[j]);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv, 0.75);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->attloc[j];
				pra1 = 0; pra2 = 0;
				if (tni->attforces[j] * vx < 0)
				{
					if (tni->attforces[j] * vz > 0)
					{
						vl = tni->attforces[j];
					}
					else {
						vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
				else
				{
					if (tni->attforces[j] * vz < 0)
					{
						vl = tni->attforces[j];
					}
					else
					{
						vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
			}
			else if (abs(tni->y_rotataion) > 20)
			{
				newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
				pra1 = 70; pra2 = 0;
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				vvx = vn^vvx;
			}

		}
		else if (s == "UL3" || s == "UR3" || s == "LL3" || s == "LR3")
		{//上下颌尖牙
			if (tni->extint > 2)
			{
				newNode->setType(AttType::optExt);//优化伸长
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (tni->extint < -2)
			{
				newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));//投影到变换后的坐标系下的x坐标
				vvx = vn^vvx;//水平矩形附件再旋转90度
				p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
				pra1 = 70; pra2 = 0;
			}
			else if (abs(tni->x_tran) > 2 || abs(tni->y_tran) > 2)
			{
				newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (abs(tni->x_rotataion) > 20 || abs(tni->rotation_theta) > 20)
			{
				newNode->setType(AttType::optRotLong);//优化旋转
				vn = NearestVer(tni, tni->attloc[j]);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->attloc[j];
				pra1 = 0; pra2 = 0;
				if (tni->attforces[j] * vx < 0)
				{
					if (tni->attforces[j] * vz > 0)
					{
						vl = tni->attforces[j];
					}
					else {
						vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
				else
				{
					if (tni->attforces[j] * vz < 0)
					{
						vl = tni->attforces[j];
					}
					else
					{
						vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
			}
			else if (abs(tni->y_rotataion) > 20)
			{
				newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
				if (tni->attforces[j] * tni->pca[0] < 0)
				{
					vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vn));
				}
				else {
					vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vn)));
				}
			}
		}
		else if (s == "UL1" || s == "UR1" || s == "UL2" || s == "UR2")//上颌切牙
		{
			if (tni->extint > 2)
			{
				newNode->setType(AttType::optExt);//优化伸长
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (tni->extint < -2)
			{
				newNode->setType(AttType::rectHori);//水平矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				ro = ArbRot(FPi / 2.0, Ou, vx)*ro;
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				vvx = vn^vvx;
				p = (tni->toothFeature.lateralPivot + tni->toothFeature.lateralCusp) / 2.0;
				pra1 = 70; pra2 = 0;
			}
			else if (abs(tni->x_tran) > 2 || abs(tni->y_tran) > 2)
			{
				newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vx.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vx.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (abs(tni->x_rotataion) > 20)
			{
				newNode->setType(AttType::optRotLong);//优化旋转
				vn = NearestVer(tni, tni->attloc[j]);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->attloc[j];
				pra1 = 0; pra2 = 0;
				if (tni->attforces[j] * vx < 0)
				{
					if (tni->attforces[j] * vz > 0)
					{
						vl = tni->attforces[j];
					}
					else {
						vl = tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
				else
				{
					if (tni->attforces[j] * vz < 0)
					{
						vl = tni->attforces[j];
					}
					else
					{
						vl = -(tni->attforces[j] - tni->attforces[j] / (1.0 / (tni->attforces[j] * vz)));
						ro = ArbRot(FPi / 3.0, Ou, vy);//绕着经过p的vv旋转rotvalue(弧度)
						vl = M(ro, vl);
					}
				}
			}
			else if (abs(tni->y_rotataion) > 20)
			{
				newNode->setType(AttType::powerRidge_b);//powerRidge
				vn=-NearestVer(tni,tni->attloc[j]);
				Point3f vv=vz.Normalize()^vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro=ArbRot(rotvalue,Ou,vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx=tni->pca[0]-vn/(1.0/(tni->pca[0]*vn));
				p=tni->attloc[j];
				pra1=0;pra2=0;
			}
		}
		else if (s == "LL1" || s == "LR1" || s == "LL2" || s == "LR2")//下颌切牙
		{
			if (abs(tni->extint) > 2)
			{
				newNode->setType(AttType::optExt);//优化伸长
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);
				Point3f vv = vz.Normalize() ^ vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro = ArbRot(rotvalue, Ou, vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
			else if (abs(tni->y_rotataion) > 20)
			{
				newNode->setType(AttType::powerRidge_b);//powerRidge
				vn=-NearestVer(tni,tni->attloc[j]);
				Point3f vv=vz.Normalize()^vn.Normalize();
				float rotvalue = acos(vz.Normalize()*vn.Normalize());
				ro=ArbRot(rotvalue,Ou,vv);//绕着经过p的vv旋转rotvalue(弧度)
				vvx=tni->pca[0]-vn/(1.0/(tni->pca[0]*vn));
				p=tni->attloc[j];
				pra1=0;pra2=0;
			}
			else
			{
				newNode->setType(AttType::rectVert);//垂直矩形附件,根据附着位置的法向量进行粘贴
				ro = ArbRot(FPi / 2.0, Ou, vy);//绕y轴旋转90度
				vn = NearestVer(tni, tni->toothFeature.lateralPivot);//找到附着面片的法向量
				Point3f vv = vx.Normalize() ^ vn.Normalize();//计算旋转轴
				float rotvalue = acos(vx.Normalize()*vn.Normalize());//计算选装角度
				ro = ArbRot(rotvalue, Ou, vv)*ro;//绕着经过p的vv旋转rotvalue(弧度)
				vvx = tni->pca[0] - vn / (1.0 / (tni->pca[0] * vn));
				p = tni->toothFeature.lateralPivot;
				pra1 = 0; pra2 = 0;
			}
		}

		if (newNode->attType<=1) continue;
		newNode->isAuto = true;
		newNode->setDuration(0,teeth_tree[current_teeth]->planstep-1);
		newNode->AttachmentNode::genAttMesh(showP,pra1,0,vl,alpha);
		if(AttType::attRootType(attype)<AttType::buttonCutout||AttType::attRootType(attype)==AttType::biteRamp||AttType::attRootType(attype)==AttType::powerPoint||AttType::attRootType(attype)==AttType::powerArm)
		{ 
			newNode->initPosition(p,vn,vvx);
		}
		newNode->updateColor(vcg::Point3f(0,0,0));
		newNode->updateRotation();
		//newNode->isAuto = true;
		newNode->vl = vl;
		newNode->vvx = vvx;
		newNode->vn = vn;
		newNode->p = p;
		teeth_tree[current_teeth]->attNodes.push_back(newNode);

		if(AttType::attRootType(newNode->attType) == AttType::optRot||AttType::attRootType(newNode->attType) == AttType::optRoot)
		{
			MeshModel *mm = md->addNewMesh(newNodeName+"0",false);
			AttachmentNode* newNNode = new AttachmentNode(mm,teeth_tree[current_teeth]->attNodes.size(),newNodeName + "0",CMN);
			newNNode->setType(newNode->attType);
			newNNode->setPickedInfo(tni->attloc[j], tni->attforces[j]);
			newNNode->setDuration(newNode->duration[0],newNode->duration[1]);
			newNNode->isAuto = true;
			newNNode->genAttMesh(showP,pra1,1,vl,alpha);
			if (AttType::attRootType(newNNode->attType)<AttType::buttonCutout||AttType::attRootType(newNNode->attType)==AttType::biteRamp||AttType::attRootType(newNNode->attType)==AttType::powerPoint)
			{
				newNNode->initPosition(p,vn,vvx);
			}
			newNNode->updateColor(vcg::Point3f(0,0,0));
			newNNode->updateRotation();
			//newNNode->isAuto = true;
			newNNode->vl = vl;
			newNNode->vvx = vvx;
			newNNode->vn = vn;
			newNNode->p = p;
			teeth_tree[current_teeth]->attNodes.push_back(newNNode);
		}
}
}