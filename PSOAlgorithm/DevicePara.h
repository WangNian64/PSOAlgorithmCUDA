#pragma once
#include <cuda_runtime.h>
#include <math.h>
#include <string>
#include <vector>
using namespace std;
//��ʾһ������
struct Vector2
{
public:
	double x;
	double y;
	//Vector2(){}
	__host__ __device__ Vector2() {}
	Vector2(double x, double y)
	{
		this->x = x;
		this->y = y;
	}
	//device���캯��
	__device__ Vector2(double x, double y, int i)//device���캯��
	{
		this->x = x;
		this->y = y;
	}
	double Distance(Vector2 v) const
	{
		return sqrt(pow(this->x - v.x, 2) + pow(this->y - v.y, 2));
	}
	bool operator<(const Vector2& v) const
	{
		if (*this == v) {
			return false;
		}
		else
		{
			if (this->x == v.x)
			{
				return this->y < v.y;
			}
			else
			{
				return this->x < v.x;
			}
		}
	}
	bool operator!=(const Vector2& v) const
	{
		return !(*this == v);
	}
	bool operator>(const Vector2& v) const
	{
		return !(*this < v || *this == v);
	}
	bool operator==(const Vector2& v) const
	{
		return abs(this->x - v.x) <= 0.0001 && abs(this->y - v.y) <= 0.0001;
	}
};
struct Vector2Int
{
public:
	int x;
	int y;
	__host__ __device__ Vector2Int() {}
	Vector2Int(int i) {}
	Vector2Int(int x, int y)
	{
		this->x = x;
		this->y = y;
	}
	__device__ Vector2Int(int x, int y, int i)
	{
		this->x = x;
		this->y = y;
	}
	__device__ double Distance(Vector2Int v) const
	{
		return sqrt(powf((this->x - v.x) / 10000.0f, 2) + powf((this->y - v.y) / 10000.0f, 2));
	}

	bool ASmallB(const Vector2Int& v) const
	{
		if (this->AEqualB(v)) {
			return false;
		}
		else
		{
			if (this->x == v.x)
			{
				return this->y < v.y;
			}
			else
			{
				return this->x < v.x;
			}
		}
	}
	bool ABigB(const Vector2Int& v) const
	{
		return !(this->ASmallB(v) || this->AEqualB(v));
	}
	bool AEqualB(const Vector2Int& v) const
	{
		return this->x == v.x && this->y == v.y;
	}
	bool ANEqualB(const Vector2Int& v) const
	{
		return !AEqualB(v);
	}


	//device�汾��
	__device__ bool ASmallB(const Vector2Int& v, int i) const
	{
		if (this->AEqualB(v, i)) {
			return false;
		}
		else
		{
			if (this->x == v.x)
			{
				return this->y < v.y;
			}
			else
			{
				return this->x < v.x;
			}
		}
	}
	__device__ bool ABigB(const Vector2Int& v, int i) const
	{
		return !(this->ASmallB(v, i) || this->AEqualB(v, i));
	}
	__device__ bool AEqualB(const Vector2Int& v, int i) const
	{
		return this->x == v.x && this->y == v.y;
	}
	__device__ bool ANEqualB(const Vector2Int& v, int i) const
	{
		return !AEqualB(v, i);
	}
};
//�豸���ࣨ��ʱ��
enum DeviceTypeEnum
{
	//�����
	HighBayWarehouse,
	DenseWarehouse,
	MultiWarehouse,

	//���ͻ�
	StraightConveryer,
	CentringConveyor,
	AlignConveyor,
	CurveConveyor,
	SpiralConveyor,
	LiftingConveyor,

	//�����������֡������е��
	LiftingTransfer,
	UniversalWheelTransfer,
	FourJointArm,

	//RGV
	RGVSys,
};
//ÿ���豸�Ķ��ز���
struct DeviceType
{
	DeviceTypeEnum typeName;//�豸������
	vector<DeviceType> adjTypes;//�ڽ��豸������
	int adjLowNum;//�豸�ڽ��豸��Ŀ�½�
	int adjUpNum;//�豸�ڽ��豸��Ŀ�Ͻ�
};
//����� ����
enum InoutType
{
	In,			//���
	Out,		//����
};

//����ڵ�ĳ���
enum PointDirect
{
	Up = 1, Right = 2, Down = 3, Left = 4
};
//�ڽӵ�ṹ
struct AdjPoint
{
public:
	int index;//�±�
	InoutType inoutType;//���������
	Vector2 pos;		//����ڵ�λ��
	PointDirect direct;//�����ķ���
	PointDirect GetDirect(string str)
	{
		PointDirect PD;
		if (str == "UP")
			PD = (PointDirect)1;
		if (str == "RIGHT")
			PD = (PointDirect)2;
		if (str == "DOWN")
			PD = (PointDirect)3;
		if (str == "LEFT")
			PD = (PointDirect)4;
		return PD;
	}
};
//�豸����˳ʱ��ת����ΪĬ��/90/180/270��
enum DeviceDirect
{
	Default,
	Rotate90,
	Rotate180,
	Rotate270
};

//�������߶εķ��򣺴�ֱ/ˮƽ
enum PathPointDirect
{
	Vert, Hori
};
//·�������Ϣ
struct PointInfo
{
	Vector2Int pointAxis;//�������
	int vertDirNum;//��Ĵ�ֱ������Ŀ
	int horiDirNum;//���ˮƽ������Ŀ
	bool isKeep;//�Ƿ���
	__device__ PointInfo() {}
	PointInfo(int i) {}
	PointInfo(Vector2Int pointAxis, int vDirNum, int hDirNum, bool iK)
	{
		this->pointAxis = pointAxis;
		vertDirNum = vDirNum;
		horiDirNum = hDirNum;
		isKeep = iK;
	}
	__device__ PointInfo(Vector2Int pointAxis, int vDirNum, int hDirNum, bool iK, int i)
	{
		this->pointAxis = pointAxis;
		vertDirNum = vDirNum;
		horiDirNum = hDirNum;
		isKeep = iK;
	}
	bool AEqualB(const PointInfo& obj)//A==B
	{
		return this->pointAxis.AEqualB(obj.pointAxis);
	}
	bool ABigB(const PointInfo& obj)//A>B
	{
		if (this->pointAxis.AEqualB(obj.pointAxis))
			return false;
		else//A!=B
		{
			if (this->pointAxis.x == obj.pointAxis.x)
			{
				return this->pointAxis.y > obj.pointAxis.y;
			}
			else
			{
				return this->pointAxis.x > obj.pointAxis.x;
			}
		}
	}
	bool ABigEqualB(const PointInfo& obj) //A>=B
	{
		return this->AEqualB(obj) || this->ABigB(obj);
	}
	bool ASmallB(const PointInfo& obj)//A<B
	{
		return !this->ABigEqualB(obj);
	}
	bool ASmallEqualB(const PointInfo& obj)//A<=B
	{
		return !this->ABigB(obj);
	}


	//device�汾�ĺ���
	__device__ bool AEqualB(const PointInfo& obj, int i)//A==B
	{
		return this->pointAxis.AEqualB(obj.pointAxis, i);
	}
	__device__ bool ABigB(const PointInfo& obj, int i)//A>B
	{
		if (this->pointAxis.AEqualB(obj.pointAxis, i))
			return false;
		else//A!=B
		{
			if (this->pointAxis.x == obj.pointAxis.x)
			{
				return this->pointAxis.y > obj.pointAxis.y;
			}
			else
			{
				return this->pointAxis.x > obj.pointAxis.x;
			}
		}
	}
	__device__ bool ABigEqualB(const PointInfo& obj, int i) //A>=B
	{
		return this->AEqualB(obj, i) || this->ABigB(obj, i);
	}
	__device__ bool ASmallB(const PointInfo& obj, int i)//A<B
	{
		return !this->ABigEqualB(obj, i);
	}
	__device__ bool ASmallEqualB(const PointInfo& obj, int i)//A<=B
	{
		return !this->ABigB(obj, i);
	}
};
//�豸�����(��0��5��������
enum DeviceRelation
{
	X, U, O, I, E, A
};
struct DeviceIDSize
{
	int ID;
	Vector2 size;
	__device__ DeviceIDSize() {}
	DeviceIDSize(int id, Vector2 s) {
		ID = id;
		size = s;
	}
	__device__ DeviceIDSize(int id, Vector2 s, int i) {
		ID = id;
		size = s;
	}
	bool operator<(const DeviceIDSize& rhs) const {
		return (this->size.x * this->size.y) < (rhs.size.x * rhs.size.y);//���������С����
	}
};
//�����豸�Ĳ���
class DevicePara
{
public:
	int ID;					//�豸ID
	double workSpeed;		//�ӹ�/����1��λ���ϵ�ʱ��
	Vector2 size;			//�豸�ߴ磨�ֱ���x���y��ĳ��ȣ�
	Vector2 axis;			//�豸����
	DeviceDirect direct;	//�豸����
	double spaceLength;		//��϶��Ϊ��ʵ�־���Լ����
	//����ڵ�����飨��Ӱ�������ߵĲ��֣�
	int adjPInCount;
	int adjPOutCount;
	AdjPoint* adjPointsIn;	//���
	AdjPoint* adjPointsOut;	//����
	DevicePara() {}
	~DevicePara() {

	}
};
//��������������ݽṹ
struct PointLink
{
public:
	int device1Index;
	int device1PointIndex;
	int device2Index;
	int device2PointIndex;
	Vector2* points;//���
	int pointNum;//�����Ŀ

	__device__ PointLink() {}
	PointLink(int i) {}
	PointLink(int device1Index, int device1PointIndex, int device2Index, int device2PointIndex, Vector2* points, int pointNum)
	{
		this->device1Index = device1Index;
		this->device1PointIndex = device1PointIndex;
		this->device2Index = device2Index;
		this->device2PointIndex = device2PointIndex;
		this->points = points;
		this->pointNum = pointNum;
	}
	__device__ PointLink(int device1Index, int device1PointIndex, int device2Index, int device2PointIndex, Vector2* points, int pointNum, int i)
	{
		this->device1Index = device1Index;
		this->device1PointIndex = device1PointIndex;
		this->device2Index = device2Index;
		this->device2PointIndex = device2PointIndex;
		this->points = points;
		this->pointNum = pointNum;
	}
};
struct InoutPoint
{
	int pointDirect;//0��ʾ��ֱ��1��ʾˮƽ
	Vector2 pointAxis;//������
	__host__ __device__ InoutPoint() {}//Ĭ�ϵ���__device__
	InoutPoint(int i) {}
	InoutPoint(int pointDirect, Vector2 pointAxis)
	{
		this->pointAxis = pointAxis;
		this->pointDirect = pointDirect;
	}
	__device__ InoutPoint(int pointDirect, Vector2 pointAxis, int i)
	{
		this->pointAxis = pointAxis;
		this->pointDirect = pointDirect;
	}
};
struct DeviceLink
{
	int inDeviceIndex;//����豸��index
	int outDeviceIndex;//�����豸��index
	int inPointIndex;//����������Ϊ�ܵ�inPoint���±꣨���ڶ���һά�ģ�
	int outPointIndex;//ͬ��
	DeviceLink() {
		inDeviceIndex = outDeviceIndex = inPointIndex = outPointIndex = 0;
	}
};
//������Ϣ
struct CargoType
{
	int deviceSum;				//�������豸��Ŀ
	int linkSum;				//�豸��Ե���Ŀ
	DeviceLink* deviceLinkList;	//�豸�����б�
	double totalVolume;			//�����ϵ���������
};
//һС��·�������ݽṹ����㣬�յ㣬���ߵķ���
struct SegPath
{
	Vector2Int p1;
	Vector2Int p2;
	PathPointDirect direct;
	__device__ SegPath() {}
	SegPath(int i) {}
	SegPath(Vector2Int p1, Vector2Int p2)
	{
		this->p1 = p1;
		this->p2 = p2;
		if (abs(p1.x - p2.x) > abs(p1.y - p2.y)) {
			direct = PathPointDirect::Hori;
		}
		else {
			direct = PathPointDirect::Vert;
		}
	}
	__device__ SegPath(Vector2Int p1, Vector2Int p2, int i)
	{
		this->p1 = p1;
		this->p2 = p2;
		if (abs(p1.x - p2.x) > abs(p1.y - p2.y)) {
			direct = PathPointDirect::Hori;
		}
		else {
			direct = PathPointDirect::Vert;
		}
	}
	bool AEqualB(const SegPath& sg)//A==B
	{
		return (this->p1.AEqualB(sg.p1) && this->p2.AEqualB(sg.p2)) || (this->p1.AEqualB(sg.p2) && this->p2.AEqualB(sg.p1));
	}
	bool ABigB(const SegPath& sg)//A>B
	{
		if ((this->p1.AEqualB(sg.p1) && this->p2.AEqualB(sg.p2)) || (this->p1.AEqualB(sg.p2) && this->p2.AEqualB(sg.p1)))
			return false;
		else//A!=B
		{
			if (this->p1.AEqualB(sg.p1))
			{
				return this->p2.ABigB(sg.p2);
			}
			else
			{
				return this->p1.ABigB(sg.p1);
			}
		}
	}
	bool ABigEqualB(const SegPath& sg) //A>=B
	{
		return this->AEqualB(sg) || this->ABigB(sg);
	}
	bool ASmallB(const SegPath& sg)//A<B
	{
		return !this->ABigEqualB(sg);
	}
	bool ASmallEqualB(const SegPath& sg)//A<=B
	{
		return !this->ABigB(sg);
	}



	/////////__device__�汾��
	__device__ bool AEqualB(const SegPath& sg, int i)//A==B
	{
		return (this->p1.AEqualB(sg.p1, i) && this->p2.AEqualB(sg.p2, i)) || (this->p1.AEqualB(sg.p2, i) && this->p2.AEqualB(sg.p1, i));
	}
	__device__ bool ABigB(const SegPath& sg, int i)//A>B
	{
		if ((this->p1.AEqualB(sg.p1, i) && this->p2.AEqualB(sg.p2, i)) || (this->p1.AEqualB(sg.p2, i) && this->p2.AEqualB(sg.p1, i)))
			return false;
		else//A!=B
		{
			if (this->p1.AEqualB(sg.p1, i))
			{
				return this->p2.ABigB(sg.p2, i);
			}
			else
			{
				return this->p1.ABigB(sg.p1, i);
			}
		}
	}
	__device__ bool ABigEqualB(const SegPath& sg, int i) //A>=B
	{
		return this->AEqualB(sg, i) || this->ABigB(sg, i);
	}
	__device__ bool ASmallB(const SegPath& sg, int i)//A<B
	{
		return !this->ABigEqualB(sg, i);
	}
	__device__ bool ASmallEqualB(const SegPath& sg, int i)//A<=B
	{
		return !this->ABigB(sg, i);
	}
};
struct StraightConveyorInfo
{
	Vector2Int startPos;
	Vector2Int endPos;
	int startVnum;
	int startHnum;
	int endVnum;
	int endHnum;
	__host__ __device__ StraightConveyorInfo() {};
	StraightConveyorInfo(int i) {};
	StraightConveyorInfo(Vector2Int sPos, Vector2Int ePos)
	{
		startPos = sPos;
		endPos = ePos;
	}
	__device__ StraightConveyorInfo(Vector2Int sPos, Vector2Int ePos, int i)
	{
		startPos = sPos;
		endPos = ePos;
	}
	bool operator<(const StraightConveyorInfo& rhs) const
	{
		if (rhs.startPos.AEqualB(this->startPos) && rhs.endPos.AEqualB(this->endPos))
		{
			return false;
		}
		else
		{
			if (this->startPos.AEqualB(rhs.startPos))
			{
				return this->endPos.ASmallB(rhs.endPos);
			}
			else
			{
				return this->startPos.ASmallB(rhs.startPos);
			}
		}
	}
};