//
//  CAstar.cpp
//  Astar
//
//  Created by xujw on 15/4/9.
//  Copyright (c) 2015年 xujw. All rights reserved.
//
//  上下左右走一格花销为10，斜着走一格花销为14，以方便计算
//  即格子宽高为10 对角线为14
#pragma once
#include "AStar.h"
//自定义排序函数
//bool mySort(const APoint* p1, const APoint* p2)
//{
//    return p1->f < p2->f;
//}


APoint::~APoint()//只存在于FitnessFunc内部
{
	if (parent) { delete[] parent; }
}

#pragma mark------CAstar-------
CAstar::~CAstar()//只存在于FitnessFunc内部
{
	if (_openList) { delete[] _openList; }
	if (_closeList) { delete[] _closeList; }
	if (_neighbourList) { delete[] _neighbourList; }
	openList_CurSize = nullptr;
	closeList_CurSize = nullptr;
	neighbourList_CurSize = nullptr;
	curPathDirect = nullptr;
}
