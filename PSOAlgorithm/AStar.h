#pragma once
//
//  CAstar.h
//  Astar
//
//  Created by xujw on 15/4/9.
//  Copyright (c) 2015年 xujw. All rights reserved.
//
/*
	F:路径评分 = g+h
	G:走一格格子的花销
	H:当前格子到目标格子的估算花销

	上下左右走一格花销为10，斜着走一格花销为14，以方便计算
	即格子宽高为10 对角线为14
 */
#ifndef __Astar__CAstar__
#define __Astar__CAstar__

 //#include "Tools.h"
#include "APoint.h"
using namespace std;

//归并排序(APoint版本)
//注意排序规则，是按照APoint的f值进行比较
static __device__ void Merge(APoint** objArray, int start, int middle, int end, APoint** tempArray)
{
	int i = start;
	int j = middle + 1;
	int index = 0;
	while (i <= middle && j <= end)
	{
		if (objArray[i]->f <= objArray[j]->f) {//排序规则
			tempArray[index++] = objArray[i++];
		}
		else {
			tempArray[index++] = objArray[j++];
		}
	}
	while (i <= middle) {
		tempArray[index++] = objArray[i++];
	}
	while (j <= end) {
		tempArray[index++] = objArray[j++];
	}
	for (int i = 0; i < index; i++) {
		objArray[start + i] = tempArray[i];
	}
}

//使用归并排序实现稳定的sort(APoint)
static __device__ void StableSort_APoint(APoint** objArray, int start, int end, APoint** tempArray)
{
	if (start < end)
	{
		int middle = (start + end) / 2;
		StableSort_APoint(objArray, start, middle, tempArray);
		StableSort_APoint(objArray, middle + 1, end, tempArray);
		Merge(objArray, start, middle, end, tempArray);
	}
}




enum PathDirection
{
	Vertical,
	Horizon
};
class CAstar
{
public:
	int openList_MaxSize;		//开放列表数组最大大小
	int closeList_MaxSize;		//关闭列表数组最大大小
	int neighbourList_MaxSize;	//周边节点列表数组最大大小
	//当前大小
	int* openList_CurSize;
	int* closeList_CurSize;
	int* neighbourList_CurSize;
	APoint** _openList;				//开放列表
	APoint** _closeList;			//关闭列表
	APoint** _neighbourList;		//周边节点
	APoint* _endPoint;
	APoint* _curPoint;

	APoint** _allPoints;			//修改这个数据结构
	int pointRowNum;				//点的行数目
	int pointColNum;				//点的列数目

	PathDirection* curPathDirect;
	__device__ CAstar() :_endPoint(nullptr), _curPoint(nullptr)//合法
	{
		openList_CurSize = new int[1];
		closeList_CurSize = new int[1];
		neighbourList_CurSize = new int[1];
		curPathDirect = new PathDirection[1];
		openList_MaxSize = 100;
		closeList_MaxSize = 100;
		neighbourList_MaxSize = 4;

		openList_CurSize[0] = 0;
		closeList_CurSize[0] = 0;
		neighbourList_CurSize[0] = 0;
	}

	~CAstar();
	//__device__ double CalcuPathLength(APoint* point);
	//__device__ APoint* findWay(PathDirection beginDirect, int beginRowIndex, int beginColIndex, int endRowIndex, int endColIndex);
	//__device__ void resetAStar();
	//bool SameDirect(APoint* curPoint, APoint* nextPoint);
private:
	//__device__ double getF(APoint* point);
	//__device__ double getH(APoint* point);
	//__device__ double getE(APoint* curPoint, APoint* nextPoint, APoint* endPoint);
	//__device__ APoint** getNeighboringPoint(APoint* point);
};

static __device__ double CalcuPathLength(APoint* point);
static __device__ APoint** getNeighboringPoint(APoint* point, int* neighbourList_CurSize, PathDirection* curPathDirect,
	int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _neighbourList);//相邻节点最多就4个
static __device__ void resetAStar(int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _openList, APoint** _closeList, APoint** _neighbourList,
	int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize);
static __device__ APoint* findWay(PathDirection* curPathDirect, PathDirection beginDirect, APoint** _allPoints, APoint* _endPoint, APoint** _neighbourList,
	APoint* _curPoint, APoint** _openList, APoint** _closeList, int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize,
	int pointColNum, int pointRowNum, int beginRowIndex, int beginColIndex, int endRowIndex, int endColIndex);
static __device__ double getF(APoint* point, APoint* _endPoint);
static __device__ double getH(APoint* point, APoint* _endPoint);
static __device__ double getE(APoint* curPoint, APoint* nextPoint, APoint* endPoint, PathDirection* curPathDirect);



//计算一条路径的总长度
static __device__ double CalcuPathLength(APoint* point)
{
	double length = 0.0;
	while (point->parent)
	{
		length += point->CalcuPointDist(*(point->parent), -1);
		point = point->parent;
	}
	return length;
}
static __device__ void resetAStar(int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _openList, APoint** _closeList, APoint** _neighbourList,
	int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize)
{
	for (int i = 0; i < pointRowNum; i++)
	{
		for (int j = 0; j < pointColNum; j++)
		{
			_allPoints[i * pointColNum + j]->resetAPoint();
		}
	}
	delete[] _openList;
	delete[] _closeList;
	delete[] _neighbourList;
	openList_CurSize[0] = 0;
	closeList_CurSize[0] = 0;
	neighbourList_CurSize[0] = 0;
}
//要改成__device__的
static __device__ APoint* findWay(PathDirection* curPathDirect, PathDirection beginDirect, APoint** _allPoints, APoint* _endPoint, APoint** _neighbourList,
	APoint* _curPoint, APoint** _openList, APoint** _closeList, int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize,
	int pointColNum, int pointRowNum, int beginRowIndex, int beginColIndex, int endRowIndex, int endColIndex)
{
	curPathDirect[0] = beginDirect;
	_endPoint = _allPoints[endRowIndex * pointColNum + endColIndex];//这个内存得是GPU内部的
	APoint* beginPoint = _allPoints[beginRowIndex * pointColNum + beginColIndex];

	if (_endPoint->type == AType::ATYPE_BARRIER)
	{
		//cout << "终点是障碍" << endl;
		return nullptr;
	}
	if (_endPoint->AEqualB(*beginPoint, -1))
	{
		_curPoint = beginPoint;
		_curPoint->parent = _endPoint;
		return _curPoint;
	}

	_openList[openList_CurSize[0]++] = beginPoint;
	beginPoint->type = AType::ATYPE_OPENED;
	beginPoint->f = getF(beginPoint, _endPoint);
	//---------
	do
	{
		//获取最小值的节点
		_curPoint = _openList[0];

		//删除最前面的节点
		if (openList_CurSize[0] == 1) {
			_openList = NULL;//这个也要改
		}
		else {
			_openList = &_openList[1];
		}
		openList_CurSize--;

		_curPoint->type = AType::ATYPE_CLOSED;
		_closeList[closeList_CurSize[0]++] = _curPoint;

		if (_curPoint->AEqualB(*_endPoint, -1))
		{
			//cout << "have find way" << endl;
			return _curPoint;
		}
		//获取相邻的节点
		APoint** neVec = getNeighboringPoint(_curPoint, neighbourList_CurSize, curPathDirect, pointRowNum, pointColNum, _allPoints, _neighbourList);
		for (int i = 0; i < neighbourList_CurSize[0]; i++)
		{
			APoint* tmpoint = neVec[i];
			double tempG = _curPoint->g + tmpoint->CalcuPointDist(*_curPoint) + getE(_curPoint, tmpoint, _endPoint, curPathDirect);
			if (tmpoint->type == AType::ATYPE_CLOSED)
			{
				continue;
			}
			//是否在开放列表里
			if (tmpoint->type != AType::ATYPE_OPENED)
			{
				tmpoint->parent = _curPoint;
				//计算GHF
				tmpoint->g = tempG;
				//更新方向
				//if (tmpoint->x == _curPoint->x)
				//{
				//    curPathDirect = PathDirection::Vertical;
				//}
				//if (tmpoint->y == _curPoint->y)
				//{
				//    curPathDirect = PathDirection::Horizon;
				//}
				tmpoint->h = getH(tmpoint, _endPoint);
				tmpoint->f = getF(tmpoint, _endPoint);
				//添加到开放列表里
				_openList[openList_CurSize[0]++] = tmpoint;
				tmpoint->type = AType::ATYPE_OPENED;
			}
			else
			{
				//已经在开放列表里
				//if (tmpoint->g < _curPoint->g)
				//if (tmpoint->h < _curPoint->h)
				if (tempG < tmpoint->g)
				{
					tmpoint->parent = _curPoint;
					tmpoint->g = tempG;
					//更新方向
					if (tmpoint->x == _curPoint->x)
					{
						curPathDirect[0] = PathDirection::Vertical;
					}
					if (tmpoint->y == _curPoint->y)
					{
						curPathDirect[0] = PathDirection::Horizon;
					}
					tmpoint->f = getF(tmpoint, _endPoint);
				}
			}
		}
		//排序 F值最小的排在前面
		//需要自己实现
		//归并排序是稳定的，用这个实现
		//stable_sort(_openList.begin(), _openList.end(), mySort);
		APoint** tempOpenList = new APoint *[openList_CurSize[0]];
		StableSort_APoint(_openList, 0, openList_CurSize[0] - 1, tempOpenList);

	} while (openList_CurSize > 0);


	//cout << "---can not find way---" << endl;
	return nullptr;
}
//bool CAstar::SameDirect(APoint* curPoint, APoint* nextPoint)
//{
//    return (curPathDirect == PathDirection::Vertical && nextPoint->x == curPoint->x)//维持原方向
//        || (curPathDirect == PathDirection::Horizon && nextPoint->y == curPoint->y);
//}
//得到F=G+H+E(E是为了对路径进行微调，减少拐点）
static __device__ double getF(APoint* point, APoint* _endPoint)
{
	return (point->g + getH(point, _endPoint));
}
//估算H
static __device__ double getH(APoint* point, APoint* _endPoint)
{
	//曼哈顿城市街区估算法
	return abs(_endPoint->y - point->y) + abs(_endPoint->x - point->x);
}
//计算E
static __device__ double getE(APoint* curPoint, APoint* nextPoint, APoint* endPoint, PathDirection* curPathDirect)
{
	//第一个点或者是直线点（不是拐点），E=0
	if (curPoint->parent == NULL)//第一个点
	{
		if (curPathDirect[0] == PathDirection::Horizon)
		{
			return (nextPoint->y == curPoint->y) ? 0 : 2;
		}
		if (curPathDirect[0] == PathDirection::Vertical)
		{
			return (nextPoint->x == curPoint->x) ? 0 : 2;
		}
	}
	if (nextPoint->x == curPoint->parent->x || nextPoint->y == curPoint->parent->y)//维持原方向
	{
		return 0.0;
	}

	//拐向终点的点
	if (nextPoint->x == endPoint->x || nextPoint->y == endPoint->y)
	{
		return 1.0;
	}
	return 2.0;
}
static __device__ APoint** getNeighboringPoint(APoint* point, int* neighbourList_CurSize, PathDirection* curPathDirect,
	int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _neighbourList)//相邻节点最多就4个
{
	neighbourList_CurSize = 0;//清空neighbor
	//可以选择根据当前方向调整点的添加顺序
	if (curPathDirect[0] == PathDirection::Vertical)//路线方向垂直，先检查垂直方向
	{
		if (point->rowIndex < pointRowNum - 1)
		{
			if (_allPoints[(point->rowIndex + 1) * pointColNum + point->colIndex]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[(point->rowIndex + 1) * pointColNum + point->colIndex];
			}
		}
		if (point->rowIndex > 0)
		{
			if (_allPoints[(point->rowIndex - 1) * pointColNum + point->colIndex]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[(point->rowIndex - 1) * pointColNum + point->colIndex];
			}
		}
		if (point->colIndex < pointColNum - 1)
		{
			if (_allPoints[point->rowIndex * pointColNum + point->colIndex + 1]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[point->rowIndex * pointColNum + point->colIndex + 1];
			}
		}
		if (point->colIndex > 0)
		{
			if (_allPoints[point->rowIndex * pointColNum + point->colIndex - 1]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[point->rowIndex * pointColNum + point->colIndex - 1];
			}
		}
	}
	else//水平方向，先检查水平方向的节点
	{
		if (point->colIndex < pointColNum - 1)
		{
			if (_allPoints[point->rowIndex * pointColNum + point->colIndex + 1]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[point->rowIndex * pointColNum + point->colIndex + 1];
			}
		}
		if (point->colIndex > 0)
		{
			if (_allPoints[point->rowIndex * pointColNum + point->colIndex - 1]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[point->rowIndex * pointColNum + point->colIndex - 1];
			}
		}
		if (point->rowIndex < pointColNum - 1)
		{
			if (_allPoints[(point->rowIndex + 1) * pointColNum + point->colIndex]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[(point->rowIndex + 1) * pointColNum + point->colIndex];
			}
		}
		if (point->rowIndex > 0)
		{
			if (_allPoints[(point->rowIndex - 1) * pointColNum + point->colIndex]->type != AType::ATYPE_BARRIER)
			{
				_neighbourList[neighbourList_CurSize[0]++] = _allPoints[(point->rowIndex - 1) * pointColNum + point->colIndex];
			}
		}
	}

	return _neighbourList;
}
#endif 