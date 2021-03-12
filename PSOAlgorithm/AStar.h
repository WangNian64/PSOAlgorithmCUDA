#pragma once
//
//  CAstar.h
//  Astar
//
//  Created by xujw on 15/4/9.
//  Copyright (c) 2015�� xujw. All rights reserved.
//
/*
	F:·������ = g+h
	G:��һ����ӵĻ���
	H:��ǰ���ӵ�Ŀ����ӵĹ��㻨��

	����������һ����Ϊ10��б����һ����Ϊ14���Է������
	�����ӿ��Ϊ10 �Խ���Ϊ14
 */
#ifndef __Astar__CAstar__
#define __Astar__CAstar__

 //#include "Tools.h"
#include "APoint.h"
using namespace std;

//�鲢����(APoint�汾)
//ע����������ǰ���APoint��fֵ���бȽ�
static __device__ void Merge(APoint** objArray, int start, int middle, int end, APoint** tempArray)
{
	int i = start;
	int j = middle + 1;
	int index = 0;
	while (i <= middle && j <= end)
	{
		if (objArray[i]->f <= objArray[j]->f) {//�������
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

//ʹ�ù鲢����ʵ���ȶ���sort(APoint)
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
	int openList_MaxSize;		//�����б���������С
	int closeList_MaxSize;		//�ر��б���������С
	int neighbourList_MaxSize;	//�ܱ߽ڵ��б���������С
	//��ǰ��С
	int* openList_CurSize;
	int* closeList_CurSize;
	int* neighbourList_CurSize;
	APoint** _openList;				//�����б�
	APoint** _closeList;			//�ر��б�
	APoint** _neighbourList;		//�ܱ߽ڵ�
	APoint* _endPoint;
	APoint* _curPoint;

	APoint** _allPoints;			//�޸�������ݽṹ
	int pointRowNum;				//�������Ŀ
	int pointColNum;				//�������Ŀ

	PathDirection* curPathDirect;
	__device__ CAstar() :_endPoint(nullptr), _curPoint(nullptr)//�Ϸ�
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
	int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _neighbourList);//���ڽڵ�����4��
static __device__ void resetAStar(int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _openList, APoint** _closeList, APoint** _neighbourList,
	int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize);
static __device__ APoint* findWay(PathDirection* curPathDirect, PathDirection beginDirect, APoint** _allPoints, APoint* _endPoint, APoint** _neighbourList,
	APoint* _curPoint, APoint** _openList, APoint** _closeList, int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize,
	int pointColNum, int pointRowNum, int beginRowIndex, int beginColIndex, int endRowIndex, int endColIndex);
static __device__ double getF(APoint* point, APoint* _endPoint);
static __device__ double getH(APoint* point, APoint* _endPoint);
static __device__ double getE(APoint* curPoint, APoint* nextPoint, APoint* endPoint, PathDirection* curPathDirect);



//����һ��·�����ܳ���
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
//Ҫ�ĳ�__device__��
static __device__ APoint* findWay(PathDirection* curPathDirect, PathDirection beginDirect, APoint** _allPoints, APoint* _endPoint, APoint** _neighbourList,
	APoint* _curPoint, APoint** _openList, APoint** _closeList, int* openList_CurSize, int* closeList_CurSize, int* neighbourList_CurSize,
	int pointColNum, int pointRowNum, int beginRowIndex, int beginColIndex, int endRowIndex, int endColIndex)
{
	curPathDirect[0] = beginDirect;
	_endPoint = _allPoints[endRowIndex * pointColNum + endColIndex];//����ڴ����GPU�ڲ���
	APoint* beginPoint = _allPoints[beginRowIndex * pointColNum + beginColIndex];

	if (_endPoint->type == AType::ATYPE_BARRIER)
	{
		//cout << "�յ����ϰ�" << endl;
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
		//��ȡ��Сֵ�Ľڵ�
		_curPoint = _openList[0];

		//ɾ����ǰ��Ľڵ�
		if (openList_CurSize[0] == 1) {
			_openList = NULL;//���ҲҪ��
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
		//��ȡ���ڵĽڵ�
		APoint** neVec = getNeighboringPoint(_curPoint, neighbourList_CurSize, curPathDirect, pointRowNum, pointColNum, _allPoints, _neighbourList);
		for (int i = 0; i < neighbourList_CurSize[0]; i++)
		{
			APoint* tmpoint = neVec[i];
			double tempG = _curPoint->g + tmpoint->CalcuPointDist(*_curPoint) + getE(_curPoint, tmpoint, _endPoint, curPathDirect);
			if (tmpoint->type == AType::ATYPE_CLOSED)
			{
				continue;
			}
			//�Ƿ��ڿ����б���
			if (tmpoint->type != AType::ATYPE_OPENED)
			{
				tmpoint->parent = _curPoint;
				//����GHF
				tmpoint->g = tempG;
				//���·���
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
				//��ӵ������б���
				_openList[openList_CurSize[0]++] = tmpoint;
				tmpoint->type = AType::ATYPE_OPENED;
			}
			else
			{
				//�Ѿ��ڿ����б���
				//if (tmpoint->g < _curPoint->g)
				//if (tmpoint->h < _curPoint->h)
				if (tempG < tmpoint->g)
				{
					tmpoint->parent = _curPoint;
					tmpoint->g = tempG;
					//���·���
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
		//���� Fֵ��С������ǰ��
		//��Ҫ�Լ�ʵ��
		//�鲢�������ȶ��ģ������ʵ��
		//stable_sort(_openList.begin(), _openList.end(), mySort);
		APoint** tempOpenList = new APoint *[openList_CurSize[0]];
		StableSort_APoint(_openList, 0, openList_CurSize[0] - 1, tempOpenList);

	} while (openList_CurSize > 0);


	//cout << "---can not find way---" << endl;
	return nullptr;
}
//bool CAstar::SameDirect(APoint* curPoint, APoint* nextPoint)
//{
//    return (curPathDirect == PathDirection::Vertical && nextPoint->x == curPoint->x)//ά��ԭ����
//        || (curPathDirect == PathDirection::Horizon && nextPoint->y == curPoint->y);
//}
//�õ�F=G+H+E(E��Ϊ�˶�·������΢�������ٹյ㣩
static __device__ double getF(APoint* point, APoint* _endPoint)
{
	return (point->g + getH(point, _endPoint));
}
//����H
static __device__ double getH(APoint* point, APoint* _endPoint)
{
	//�����ٳ��н������㷨
	return abs(_endPoint->y - point->y) + abs(_endPoint->x - point->x);
}
//����E
static __device__ double getE(APoint* curPoint, APoint* nextPoint, APoint* endPoint, PathDirection* curPathDirect)
{
	//��һ���������ֱ�ߵ㣨���ǹյ㣩��E=0
	if (curPoint->parent == NULL)//��һ����
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
	if (nextPoint->x == curPoint->parent->x || nextPoint->y == curPoint->parent->y)//ά��ԭ����
	{
		return 0.0;
	}

	//�����յ�ĵ�
	if (nextPoint->x == endPoint->x || nextPoint->y == endPoint->y)
	{
		return 1.0;
	}
	return 2.0;
}
static __device__ APoint** getNeighboringPoint(APoint* point, int* neighbourList_CurSize, PathDirection* curPathDirect,
	int pointRowNum, int pointColNum, APoint** _allPoints, APoint** _neighbourList)//���ڽڵ�����4��
{
	neighbourList_CurSize = 0;//���neighbor
	//����ѡ����ݵ�ǰ�������������˳��
	if (curPathDirect[0] == PathDirection::Vertical)//·�߷���ֱ���ȼ�鴹ֱ����
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
	else//ˮƽ�����ȼ��ˮƽ����Ľڵ�
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