//
//  CAstar.cpp
//  Astar
//
//  Created by xujw on 15/4/9.
//  Copyright (c) 2015�� xujw. All rights reserved.
//
//  ����������һ����Ϊ10��б����һ����Ϊ14���Է������
//  �����ӿ��Ϊ10 �Խ���Ϊ14
#pragma once
#include "AStar.h"
//�Զ���������
//bool mySort(const APoint* p1, const APoint* p2)
//{
//    return p1->f < p2->f;
//}


APoint::~APoint()//ֻ������FitnessFunc�ڲ�
{
	if (parent) { delete[] parent; }
}

#pragma mark------CAstar-------
CAstar::~CAstar()//ֻ������FitnessFunc�ڲ�
{
	if (_openList) { delete[] _openList; }
	if (_closeList) { delete[] _closeList; }
	if (_neighbourList) { delete[] _neighbourList; }
	openList_CurSize = nullptr;
	closeList_CurSize = nullptr;
	neighbourList_CurSize = nullptr;
	curPathDirect = nullptr;
}
