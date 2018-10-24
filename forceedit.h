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
#ifndef FORCEEDITPLUGIN_H
#define FORCEEDITPLUGIN_H

#include <QObject>
#include <common/interfaces.h>
#include <common/AttNameDefine.h>
#include <common/Attachmodel.h>
#include <common/Mq_teethtree.h>
#include <meshlab/glarea.h>
#include <wrap/gui/rubberband.h>
#include "forcewidget.h"

//class minCFT
//{
//public:
//	vcg::Point3f bp;
//	vcg::Point3f Fa;
//	vcg::Point3f Ta;
//};

namespace vcg{
	enum ForceMark{
		FMY,FMN
	};
};

struct ForceSet {
	vcg::Point3f f_beginpos;
	vcg::Point3f f_direction;
	double f_value;
};


class ForceEditPlugin : public QObject, public MeshEditInterface
{
	Q_OBJECT
	Q_INTERFACES(MeshEditInterface)
		
public:

	static void Mark(CVertexO *v,vcg::ForceMark m){
		v->IMark() = m;
	}

	ForceEditPlugin();
	virtual ~ForceEditPlugin(){}
	static const QString Info();
	bool StartEdit(MeshDocument& _md,GLArea*,TeethTree a[],bool has_prepare);
	void EndEdit(MeshModel& /*m*/,GLArea* /*parent*/);
	void Decorate(MeshModel& /*m*/,GLArea* /*parent*/,QPainter* p);
	void Decorate(MeshModel& /*m*/,GLArea*){};
	void mousePressEvent(QMouseEvent*,MeshModel&,GLArea*);
	void mouseMoveEvent(QMouseEvent*,MeshModel&,GLArea*);
	void mouseReleaseEvent(QMouseEvent* event,MeshModel& /*m*/,GLArea*);
	virtual void LayerChanged(MeshDocument& md,MeshModel& oldMeshModel,GLArea* parent){};

public:
	void getCurrentManuMode();
	void DrawManipulators(ToothNode* node,GLArea* gla);
	void DrawTranslateManipulators(ToothNode* node,GLArea* gla);
	void DrawRotateManipulators(ToothNode* node,GLArea *gla);

	void ChangeCurrentTeethTree(int index){
		this->current_teeth = index;
		CMN = teeth_tree[current_teeth]->nodeList[0];
		foreach(ToothNode *mni,teeth_tree[current_teeth]->nodeList)
		{
			mni->m->visible=true;
			emit updatemeshvisiable();
		}
	};
	ToothNode* CurrentNode() {return CMN;}
	void SetCurrentNode(ToothNode* cn)
	{
		CMN = cn;
		if (CMN->isattachment)
		{
			forceui->ui.lb_attName->setText(CMN->qnodename);
		}
		else
			forceui->ui.lb_attName->setText("XXXX");
	};
	void deleteAtt(AttachmentNode* node);

public:
	vcg::Matrix44f original_Transform;
	vcg::Matrix44f delta_Transform;



	/////// 2018.5.25 起 - 计算所有期的带附件牙齿以及牙套，存入vector中 - GY
	//void slotAttachOptim();
	std::vector<vector<ToothNode*>> all_Braces;
	bool twoTeethCollisionDetecte(ToothNode* td1, ToothNode* td2);	// 两颗牙齿检测碰撞
	void twoTeethCollisionRemove(ToothNode* td1, ToothNode* td2);	// 碰撞处理
	int compute_all_Teeth_and_Braces(int start, int end, int stride);
	void compute_all_Force_and_Torque(std::map<std::string, std::vector<float>>& error_map, int end);
	void ComputeAlignerForceofTooth(std::map<CVertexO *,double> &res);
	//int max_step; // 记录上下颌牙齿治疗期数的最大值
	/////// 2018.5.25 起 - 计算所有期的带附件牙齿以及牙套，存入vector中 - GY



    int current_teeth; //0：max 1:men
    std::vector<TeethTree *> teeth_tree;
	MeshDocument* md;
	GLArea *gla;
    enum ManuMode
	{
		TranX,
		TranY,
		TranZ,
        RotX,	
		RotY,
		RotZ,
		ScaleX,
		ScaleY,
		ScaleZ,
		None
    };
	std::vector<vcg::Point3f> showP;
	
	void slotComputeForceofsteps(int i, int current_teeth);
	void ComputeAttachForOneTooth(ToothNode* tni);
private:
	//保存牙套与牙齿之间的力
	vcg::Point3f ForceVector;

	bool forceJudge;
    bool has_init;
	bool pickloc;
	bool isRightClicked;
	bool hasnewfata;

	bool is_measure,was_ready;
	bool firstclick;
	vcg::Rubberband rubberband;
	
    ForceWidget* forceui;
    ToothNode* CMN;

	int attype;
	ManuMode curMode;
	int importID;

	double inputValue;
	float AxisSize;
	QPoint curpos;//记录放置附件的位置
	QPoint curRightPos;

	//附件的移动与旋转
	bool isMoving;
	vcg::Point2i startdrag;
	vcg::Point2i enddrag;

	float currScreenOffset_X;
	float currScreenOffset_Y;

	float displayOffset;
	float displayOffset_X;
	float displayOffset_Y;
	float displayOffset_Z;

	float currOffset;
	float currOffset_X;
	float currOffset_Y;
	float currOffset_Z;

	vcg::Point3f screen_xaxis;
	vcg::Point3f screen_yaxis;
	vcg::Point3f screen_zaxis;

private:
	void DrawWireCone(vcg::Point3f loc,vcg::Point3f direct);
	void DrawVector(vcg::Point3f &vector,vcg::Color4b color);
	void DrawCircle(float r, float g, float b);
	void DrawArrows(float r, float g, float b);
	void DrawCubes(float r, float g, float b);
	void UpdateMatrix(MeshModel& model, GLArea * gla, bool applymouseoffset, bool useinputnumber=false);
	void DrawMeshBox(MeshModel& model);
	void BuildAttach(CFaceO* face);

	void resetOffsets();
	void applyMotion(MeshModel &model, GLArea *gla);
	void cancelMotion(MeshModel &model, GLArea *gla);
	void FreezeTransfrom();

	bool MyPick(const int &x, const int &y, vcg::Point3f &pp, float mydepth);
	void popAttAdjMenu(QPoint pos);

	void Calattloc(bool,vcg::Point3f);

	void importAttFile();

	vcg::Point3f NearestVer(ToothNode *tni,vcg::Point3f v);
	vcg::Point3f NearestVer(ToothNode *tni,vcg::Point3f v,vcg::Point3f xa);
	vcg::Matrix44f ArbRot(float t,vcg::Point3f p, vcg::Point3f vv,float alpha=1);

public slots:
	void slotImportAttInfo();
	void slotExportAttInfo();
	void slotExportModel();
	void slotCalRefTorque();
	void slotMeasure();
	void slotAutoAtt();
	void slotBuildAtt(int);
	void slotGetManuValue(double);
	void slotManStart();
	void slotManApply();
	void slotChangeDuration();
	void slotChangeAttType();
	void slotAttachOptim();
	void slotAdjustParam(int);
	void slotCalDesiredForces();
	//旧版本的导入和导出
	//void slotImportOld();
	//void slotExportOld();

signals:
	void suspendEditToggle();
	void updatemeshvisiable();

};

#endif
