#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
enum AType
{
	ATYPE_UNKNOWN,
	ATYPE_CLOSED,
	ATYPE_OPENED,
	ATYPE_BARRIER   //障碍
};

class APoint
{
public:
	__device__ APoint() :x(0)
		, y(0)
		, h(0)
		, f(0)
		, g(0)
		, parent(nullptr)
		, type(AType::ATYPE_UNKNOWN)
	{

	}
	~APoint();
	int rowIndex;
	int colIndex;
	double x;
	double y;
	AType type;		//类型:障碍、开放列表、关闭列表
	double f;		//f = g+h
	double g;
	double h;
	APoint* parent;
	bool AEqualB(const APoint& po)
	{
		if (x == po.x && y == po.y)
		{
			return true;
		}
		return false;
	}
	__device__ bool AEqualB(const APoint& po, int i)
	{
		if (x == po.x && y == po.y)
		{
			return true;
		}
		return false;
	}
	__device__ double CalcuPointDist(const APoint& point)
	{
		return abs(this->y - point.y) + abs(this->x - point.x);
	}
	__device__ double CalcuPointDist(const APoint& point, int i)
	{
		return abs(this->y - point.y) + abs(this->x - point.x);
	}
	__device__ void resetAPoint()
	{
		f = g = h = 0;
		type = AType::ATYPE_UNKNOWN;
		parent = nullptr;
	}
};