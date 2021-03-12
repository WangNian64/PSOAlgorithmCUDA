#pragma once
#include <algorithm>
#include <math.h> 
#include "Tools.h"
#include "AStar.h"
#include "PSO.h"
#pragma region �ж�������������»������ҹ�ϵ
static __device__ bool IsInLeft2Out(Vector2 inPos, Vector2 outPos)
{
	return inPos.x <= outPos.x;
}
static __device__ bool IsInRight2Out(Vector2 inPos, Vector2 outPos)
{
	return inPos.x >= outPos.x;
}
static __device__ bool IsInUp2Out(Vector2 inPos, Vector2 outPos)
{
	return inPos.y >= outPos.y;
}
static __device__ bool IsInDown2Out(Vector2 inPos, Vector2 outPos)
{
	return inPos.y <= outPos.y;
}

#pragma endregion
static __global__ void FitnessFunction(int curIterNum, int maxIterNum, int particleNum, int* bestParticleIndex,
	/*ProblemParas proParas, �̶������ģ����ù�*/
	int DeviceSum, int fixedLinkPointSum, int fixedUniqueLinkPointSum, int vertPointCount, int horiPointCount, double workShopLength, double workShopWidth, double convey2DeviceDist, /*double conveyWidth, */
	double strConveyorUnitCost, double curveConveyorUnitCost, double conveyMinDist, /*double conveyMinLength, */double conveySpeed, Vector2 entrancePos, Vector2 exitPos,
	int CargoTypeNum, int totalLinkSum,

	/*CargoType* �̶�����*/
	/*int* deviceSum, */int* linkSum, int* accumLinkSum, DeviceLink* deviceLinkList, double* totalVolume,

	/*DevicePara**/
	Vector2* size, double* spaceLength, int* adjPInCount, int* adjPOutCount, int* accumAdjPInCount, int* accumAdjPOutCount,
	int totalInPoint, int totalOutPoint, AdjPoint* adjPointsIn, AdjPoint* adjPointsOut,
	/*Particle*/
	int dim, int fitnessCount, double* fitness_GPU, double* position_GPU, /*double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU*/
	/*�洢��������������·��Ϣ*/
	double* curBestFitnessVal, int inoutPSize, InoutPoint* inoutPoints, StraightConveyorInfo* strConveyorList,
	int* strConveyorListSum, Vector2Int* curveConveyorList, int* curveConveyorListSum,
	int* pointDirectArray, curandState* globalState);

double CalcuTotalArea(Particle& particle, DevicePara* copyDeviceParas);
//double CalcuDeviceDist(Vector2 pos1, Vector2 pos2);
static __device__ int FindAxisIndex(double axis, const double* axisList, int axisCount);

//˳ʱ����ת�������
static __device__ Vector2 Rotate(Vector2 pointPos, Vector2 centerPos, float rotateAngle);
//��һ������*10000Ȼ���������뵽int
static __device__ int Multi10000ToInt(double num);

static __device__ Vector2Int Multi10000ToInt(Vector2 v);
//��Ӧ�ȼ��㺯�� GPU
static __global__ void FitnessFunction(int curIterNum, int maxIterNum, int particleNum, int* bestParticleIndex,
	/*ProblemParas proParas, �̶������ģ����ù�*/
	int DeviceSum, int fixedLinkPointSum, int fixedUniqueLinkPointSum, int vertPointCount, int horiPointCount, double workShopLength, double workShopWidth, double convey2DeviceDist, /*double conveyWidth, */
	double strConveyorUnitCost, double curveConveyorUnitCost, double conveyMinDist, /*double conveyMinLength, */double conveySpeed, Vector2 entrancePos, Vector2 exitPos,
	int CargoTypeNum, int totalLinkSum,

	/*CargoType* �̶�����*/
	/*int* deviceSum, */int* linkSum, int* accumLinkSum, DeviceLink* deviceLinkList, double* totalVolume,

	/*DevicePara**/
	Vector2* size, double* spaceLength, int* adjPInCount, int* adjPOutCount, int* accumAdjPInCount, int* accumAdjPOutCount,
	int totalInPoint, int totalOutPoint, AdjPoint* adjPointsIn, AdjPoint* adjPointsOut,
	/*Particle*/
	int dim, int fitnessCount, double* fitness_GPU, double* position_GPU, /*double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU*/
	/*�洢��������������·��Ϣ*/
	double* curBestFitnessVal, int inoutPSize, InoutPoint* inoutPoints, StraightConveyorInfo* strConveyorList,
	int* strConveyorListSum, Vector2Int* curveConveyorList, int* curveConveyorListSum,
	int* pointDirectArray, curandState* globalState)
{
	//���ӵ��±�i��Ҫ�Լ�����
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	double punishValue1 = 0;
	double punishValue2 = 0;
	bool IsDeviceOverlap = false;//�Ƿ��ص�
	double deviceDist = 0;
	fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = 0;

#pragma region ���һ���豸����
	//��Ŀ��DeviceSum
	//DevicePara* copyDeviceParas_copy = new DevicePara[DeviceSum];
	//ֻ��size��Ҫcopy
	Vector2* size_copy = new Vector2[DeviceSum];					//�豸�ߴ磨�ֱ���x���y��ĳ��ȣ�

	//ֻ��size��Ҫ����һ���µ�
	for (int i = 0; i < DeviceSum; i++)
	{
		size_copy[i] = size[i];
	}
#pragma endregion

#pragma region �����豸���򣬵����豸�ߴ�xy�ͳ��������
	int curAdjPIn_Index = 0;
	int curAdjPOut_Index = 0;
	for (int i = 2; i < dim; i += 3)
	{
		//doubleתint��ת��ΪDirection��Ȼ����ݳ������¼����豸�ߴ�ͳ����
		//Rotate90����Rotate270���ߴ��x��y����
		//����ڰ���˳ʱ����,��ת��=Direction*90(���ö�Ӧ0,90,180,270��
		DeviceDirect curDirect = (DeviceDirect)(int)position_GPU[index * dim + i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			double tempX = size_copy[i / 3].x;
			size_copy[i / 3].x = size_copy[i / 3].y;
			size_copy[i / 3].y = tempX;
		}
		//���¼�����ת��ĳ��������
		Vector2 deviceCenterPos(0, 0, -1);//��֤���õ���device�汾��
		double rotateAngle = curDirect * 90;
		int newDirect = 0;

		for (int pointIndex = 0; pointIndex < adjPInCount[i / 3]; ++pointIndex)
		{
			AdjPoint& point = adjPointsIn[curAdjPIn_Index++];
			point.pos = Rotate(point.pos, deviceCenterPos, rotateAngle);
			newDirect = point.direct + (int)curDirect;
			point.direct = (newDirect == 4) ? (PointDirect)4 : (PointDirect)(newDirect % 4);
		}
		for (int pointIndex = 0; pointIndex < adjPOutCount[i / 3]; ++pointIndex)
		{
			AdjPoint& point = adjPointsOut[curAdjPOut_Index++];
			point.pos = Rotate(point.pos, deviceCenterPos, rotateAngle);
			newDirect = point.direct + (int)curDirect;
			point.direct = (newDirect == 4) ? (PointDirect)4 : (PointDirect)(newDirect % 4);
		}
	}
#pragma endregion

#pragma region ����豸�Ƿ��ص�
	//����ص������е���
	//���ͱ�׼�ᷢ��ʲô��
	double outSizeLength, outSizeWidth;
	for (int i = 0; i < dim; i += 3) {
		outSizeLength = 0.5 * size_copy[i / 3].x + spaceLength[i / 3];
		outSizeWidth = 0.5 * size_copy[i / 3].y + spaceLength[i / 3];
		double firstLowX = position_GPU[index * dim + i] - outSizeLength;
		double firstUpX = position_GPU[index * dim + i] + outSizeLength;
		double firstLowY = position_GPU[index * dim + i + 1] - outSizeWidth;
		double firstUpY = position_GPU[index * dim + i + 1] + outSizeWidth;
		for (int j = i + 3; j < dim; j += 3) {
			outSizeLength = 0.5 * size_copy[j / 3].x + spaceLength[j / 3];
			outSizeWidth = 0.5 * size_copy[j / 3].y + spaceLength[j / 3];
			double secondLowX = position_GPU[index * dim + j] - outSizeLength;
			double secondUpX = position_GPU[index * dim + j] + outSizeLength;
			double secondLowY = position_GPU[index * dim + j + 1] - outSizeWidth;
			double secondUpY = position_GPU[index * dim + j + 1] + outSizeWidth;
			if (IsRangeOverlap(firstLowX, firstUpX, secondLowX, secondUpX) && IsRangeOverlap(firstLowY, firstUpY, secondLowY, secondUpY)) {
				//particle.fitness_[0] = particle.fitness_[1] = MAX_FITNESS;
				IsDeviceOverlap = true;
				//cout << curIterNum << ":" << "�ص�" << endl;
				//return;
			}
		}
	}
#pragma endregion

#pragma region ����豸�ص��������豸λ�ã����������ص�������ֱ�Ӳ���

#pragma region �����豸λ��
	if (IsDeviceOverlap == true) {
		//1.��һ���豸�ĳߴ�����洢�����豸
		//2.ÿ�δӵ�һ���豸��ʼ�����Ҵӱ��Ϊ1���豸��ʼ����Ƿ�͸��豸�ص�
		//3.����ص���ִ�����²�����
		//	a.�涨�豸ֻ���������ƶ�
		//  b.�ֱ������豸����/���ƶ��ľ���
		//	c.ԭ��ѡ���ƶ��������豸�ĳߴ�x֮��/�ߴ�y֮��������ڳ���xy�ߴ������С�ķ����ƶ�
		//Ŀǰ���ƶ���������λ�ã��Ȳ��������
		//4.�ƶ�֮���´λ��Ǵӵ�һ���豸��ʼ��⣬ֱ�����еĶ����ص�
		//5.Ȼ�����������豸λ�õİ�����Σ����С�ڳ���ߴ磬�ƶ��������ڵ�һ�����λ��
		//Ȼ���������ƶ�ǰ���λ�ò��޸������豸��λ��
		//�ñ�����������������ر��ף���������ά��һ����İ������
		int deviceIDSizeCount = DeviceSum;
		DeviceIDSize* deviceIDSizeList = new DeviceIDSize[DeviceSum];//�����豸��С�����ID����
		//���������ṹ���豸���꣬��Ϊ���ܻ��޸�ʧ��
		double* particlePosList = new double[dim];
		for (int i = 0; i < dim; ++i) {
			particlePosList[i] = position_GPU[index * dim + i];
		}
		for (int i = 0; i < deviceIDSizeCount; ++i) {
			deviceIDSizeList[i] = DeviceIDSize(i, size_copy[i], -1);
		}
		DeviceIDSize_Sort(deviceIDSizeList, 0, deviceIDSizeCount - 1);//�����豸�ĳߴ�����

		double outSizeLength1, outSizeWidth1;
		double outSizeLength2, outSizeWidth2;
		int firstID, secondID;
		int maxIter = 1000;//��ֹ��ѭ��
		int curIter = 0;
		bool tooMuch = false;
		for (int i = 0; i < deviceIDSizeCount; ++i) {
			//����������豸�Ƿ�����ص�
			firstID = deviceIDSizeList[i].ID;
			outSizeLength1 = 0.5 * size_copy[firstID].x + spaceLength[firstID];
			outSizeWidth1 = 0.5 * size_copy[firstID].y + spaceLength[firstID];
			double firstLowX = particlePosList[3 * firstID] - outSizeLength1;
			double firstUpX = particlePosList[3 * firstID] + outSizeLength1;
			double firstLowY = particlePosList[3 * firstID + 1] - outSizeWidth1;
			double firstUpY = particlePosList[3 * firstID + 1] + outSizeWidth1;
			for (int j = 0; j < deviceIDSizeCount;) {
				++curIter;
				if (curIter > maxIter) {
					fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
					printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
					return;
				}
				secondID = deviceIDSizeList[j].ID;
				if (firstID != secondID) {
					outSizeLength2 = 0.5 * size_copy[secondID].x + spaceLength[secondID];
					outSizeWidth2 = 0.5 * size_copy[secondID].y + spaceLength[secondID];
					double secondLowX = particlePosList[3 * secondID] - outSizeLength2;
					double secondUpX = particlePosList[3 * secondID] + outSizeLength2;
					double secondLowY = particlePosList[3 * secondID + 1] - outSizeWidth2;
					double secondUpY = particlePosList[3 * secondID + 1] + outSizeWidth2;
					if (IsRangeOverlap(firstLowX, firstUpX, secondLowX, secondUpX)
						&& IsRangeOverlap(firstLowY, firstUpY, secondLowY, secondUpY)) {
						j = 0;//ֻҪ�����ص��ľ�Ҫ���¿�ʼ�ж�
						//�ص��ˣ��ƶ��豸�����ݱ���ȷ�������ƶ�
						//������Բ�һ���ܺã���Ҫ��֤
						double rateLeft = (outSizeLength1 + outSizeLength2) * 2 / workShopLength;
						double rateDown = (outSizeWidth1 + outSizeWidth2) * 2 / workShopWidth;
						if (rateLeft < rateDown) {//˵��Ӧ�������ƶ�
							particlePosList[firstID * 3] = secondLowX - outSizeLength1 - 0.1;
						}
						else {//�����ƶ�
							particlePosList[firstID * 3 + 1] = secondLowY - outSizeWidth1 - 0.1;
						}
						firstLowX = particlePosList[3 * firstID] - outSizeLength1;
						firstUpX = particlePosList[3 * firstID] + outSizeLength1;
						firstLowY = particlePosList[3 * firstID + 1] - outSizeWidth1;
						firstUpY = particlePosList[3 * firstID + 1] + outSizeWidth1;
					}
					else {
						++j;
					}
				}
				else {
					++j;
				}
			}
		}
		//5.Ȼ�����������豸λ�õİ�����Σ����С�ڳ���ߴ磬�ƶ��������ڵ�һ�����λ��
		//ʵ�֣�����������豸���ĸ�����ı߽�
		double min_X, min_Y, max_X, max_Y;
		min_X = min_Y = INT_MAX;
		max_X = max_Y = -INT_MAX;
		for (int i = 0; i < DeviceSum; ++i) {
			double outSizeLength = size_copy[i].x * 0.5 + spaceLength[i];
			double outSizeWidth = size_copy[i].y * 0.5 + spaceLength[i];
			min_X = getMin(min_X, particlePosList[3 * i] - outSizeLength);
			max_X = getMax(max_X, particlePosList[3 * i] + outSizeLength);
			min_Y = getMin(min_Y, particlePosList[3 * i + 1] - outSizeWidth);
			max_Y = getMax(max_Y, particlePosList[3 * i + 1] + outSizeWidth);
		}
		Vector2 oriRectAxis((max_X + min_X) / 2.0, (max_Y + min_Y) / 2.0, -1);//�ܰ�����ε���������
		Vector2 newRectAxis;
		if ((max_X - min_X) <= workShopLength
			&& (max_Y - min_Y) <= workShopWidth) {//�������С�ڳ����С
			//��������������һ������
			double rectLowX = 0 + (max_X - min_X) * 0.5;
			double rectHighX = workShopLength - (max_X - min_X) * 0.5;
			double rectLowY = 0 + (max_Y - min_Y) * 0.5;
			double rectHighY = workShopWidth - (max_Y - min_Y) * 0.5;
			newRectAxis.x = createARandomNum(globalState, index) * (rectHighX - rectLowX) + rectLowX;
			newRectAxis.y = createARandomNum(globalState, index) * (rectHighY - rectLowY) + rectLowY;
			//���ݰ�����ε�����仯���޸������豸������
			double deviceOffsetX = newRectAxis.x - oriRectAxis.x;
			double deviceOffsetY = newRectAxis.y - oriRectAxis.y;
			for (int i = 0; i < DeviceSum; ++i) {
				particlePosList[3 * i] += deviceOffsetX;
				particlePosList[3 * i + 1] += deviceOffsetY;
			}
			IsDeviceOverlap = false;
			//�޸�postion����Ҫ�޸���������
			for (int i = 0; i < DeviceSum; i++) {
				position_GPU[index * dim + 3 * i] = particlePosList[3 * i];
				position_GPU[index * dim + 3 * i + 1] = particlePosList[3 * i + 1];
			}
		}
		else {
			IsDeviceOverlap = true;
			fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
			printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
			return;
		}
	}
#pragma endregion

	if (IsDeviceOverlap == false)
	{
		printf("ewfweff");
#pragma region ���뷽��1�������豸���ĵ�x��y
		////if (curIterNum == 199)
		////{
		//for (int i = 0; i < proParas.DeviceSum; i++)
		//{
		//	for (int j = 0; j < proParas.DeviceSum; j++)
		//	{
		//		if (i != j)
		//		{
		//			if (abs(particle.position_[3 * i] - particle.position_[3 * j]) <= 1)
		//			{
		//				particle.position_[3 * i] = particle.position_[3 * j]
		//					= (particle.position_[3 * i] + particle.position_[3 * j]) * 0.5f;
		//			}
		//			if (abs(particle.position_[3 * i + 1] - particle.position_[3 * j + 1]) <= 1)
		//			{
		//				particle.position_[3 * i + 1] = particle.position_[3 * j + 1]
		//					= (particle.position_[3 * i + 1] + particle.position_[3 * j + 1]) * 0.5f;
		//			}
		//		}
		//	}
		//}
		////�������
		//for (int i = 0; i < proParas.DeviceSum; i++)
		//{
		//	if (abs(particle.position_[3 * i] - proParas.entrancePos.x) <= 1)
		//	{
		//		particle.position_[3 * i] = proParas.entrancePos.x;
		//	}
		//	if (abs(particle.position_[3 * i + 1] - proParas.entrancePos.y) <= 1)
		//	{
		//		particle.position_[3 * i + 1] = proParas.entrancePos.y;
		//	}
		//}
		////}
#pragma endregion

#pragma region ���뷽��2������豸����ڵ����겢���ж������
		//��������cargoTypeList
		int curDeviceLink_Index = 0;//�ܵ��±�
		for (int j = 0; j < CargoTypeNum; j++)
		{
			for (int k = 0; k < linkSum[j]; k++)
			{
				int outDeviceIndex = deviceLinkList[curDeviceLink_Index].outDeviceIndex;
				int inDeviceIndex = deviceLinkList[curDeviceLink_Index].inDeviceIndex;
				int outPointIndex = deviceLinkList[curDeviceLink_Index].outPointIndex;
				int inPointIndex = deviceLinkList[curDeviceLink_Index].inPointIndex;
				AdjPoint outPoint, inPoint;
				//�������1���ֿ����
				if (outDeviceIndex == -1)
				{
					//�����и����⣺�±�Ķ�λ
					//ÿ��adjPoint����Ĵ�С��һ���������������һ����¼ǰ��x��adj�������Ŀ������
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (inPointTPos.x != entrancePos.x && abs(inPointTPos.x - entrancePos.x) < conveyMinDist)
					{
						//ֻ���޸�in�������޸����
						double moveLength = inPointTPos.x - entrancePos.x;
						position_GPU[index * dim + inDeviceIndex * 3] -= moveLength;

					}
					else if (inPointTPos.y != entrancePos.y && abs(inPointTPos.y - entrancePos.y) < conveyMinDist)
					{
						double moveLength = inPointTPos.y - entrancePos.y;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else if (inDeviceIndex == -2)//�������2���ֿ����
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != exitPos.x && abs(outPointTPos.x - exitPos.x) < conveyMinDist)
					{
						//ֻ���޸�out�������޸ĳ���
						double moveLength = outPointTPos.x - exitPos.x;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
					}
					else if (outPointTPos.y != exitPos.y && abs(outPointTPos.y - exitPos.y) < conveyMinDist)
					{
						double moveLength = outPointTPos.y - exitPos.y;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else//�������
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					//�ڿ������豸���������¶Ա�
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
					{
						//x����ӽ�
						double moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
						position_GPU[index * dim + inDeviceIndex * 3] += moveLength;

					}
					else if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)
					{
						//y����ӽ�
						double moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] += moveLength;
					}
				}

				curDeviceLink_Index++;
			}
		}
#pragma endregion

#pragma region ���뷽��3��������Եĳ���ڵ�ĳ��������Ż�
		//��������cargoTypeList
		//if ��һ���ǳ���Ϊi�豸���Ҳ������һ������ô�Ϳ����ó���һ�Գ���ڵ�
		curDeviceLink_Index = 0;
		for (int j = 0; j < CargoTypeNum; j++)
		{
			for (int k = 0; k < linkSum[j]; k++)
			{
				int outDeviceIndex = deviceLinkList[curDeviceLink_Index].outDeviceIndex;
				int inDeviceIndex = deviceLinkList[curDeviceLink_Index].inDeviceIndex;
				int outPointIndex = deviceLinkList[curDeviceLink_Index].outPointIndex;
				int inPointIndex = deviceLinkList[curDeviceLink_Index].inPointIndex;
				AdjPoint outPoint, inPoint;
				PointDirect outPointDirect, inPointDirect;
				//�������1���ֿ����
				if (outDeviceIndex == -1)
				{
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (inPointTPos.x != entrancePos.x && abs(inPointTPos.x - entrancePos.x) < conveyMinDist)
					{
						//ֻ���޸�in�������޸����
						double moveLength = inPointTPos.x - entrancePos.x;
						position_GPU[index * dim + inDeviceIndex * 3] -= moveLength;

					}
					else if (inPointTPos.y != entrancePos.y && abs(inPointTPos.y - entrancePos.y) < conveyMinDist)
					{
						double moveLength = inPointTPos.y - entrancePos.y;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else if (inDeviceIndex == -2)//�������2���ֿ����
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != exitPos.x && abs(outPointTPos.x - exitPos.x) < conveyMinDist)
					{
						//ֻ���޸�out�������޸ĳ���
						double moveLength = outPointTPos.x - exitPos.x;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
					}
					else if (outPointTPos.y != exitPos.y && abs(outPointTPos.y - exitPos.y) < conveyMinDist)
					{
						double moveLength = outPointTPos.y - exitPos.y;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else//�������
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					//���ǵ�ĳ���
					outPointDirect = outPoint.direct;
					inPointDirect = inPoint.direct;
					//����ڵ���ʵ����
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					//������豸�ĸ����򲻿������ı߽�
					Vector2 outDeviceUpPos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1] + 0.5 * size_copy[outDeviceIndex].y, -1);
					Vector2 outDeviceDownPos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1] - 0.5 * size_copy[outDeviceIndex].y, -1);
					Vector2 outDeviceLeftPos(position_GPU[index * dim + outDeviceIndex * 3] - 0.5 * size_copy[outDeviceIndex].x, position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 outDeviceRightPos(position_GPU[index * dim + outDeviceIndex * 3] + 0.5 * size_copy[outDeviceIndex].x, position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);

					Vector2 inDeviceUpPos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1] + 0.5 * size_copy[inDeviceIndex].y, -1);
					Vector2 inDeviceDownPos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1] - 0.5 * size_copy[inDeviceIndex].y, -1);
					Vector2 inDeviceLeftPos(position_GPU[index * dim + inDeviceIndex * 3] - 0.5 * size_copy[inDeviceIndex].x, position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					Vector2 inDeviceRightPos(position_GPU[index * dim + inDeviceIndex * 3] + 0.5 * size_copy[inDeviceIndex].x, position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);


					//������豸�ĸ�����ı߽�
					Vector2 outUpPos(outDeviceUpPos.x, outDeviceUpPos.y + convey2DeviceDist, -1);
					Vector2 outDownPos(outDeviceDownPos.x, outDeviceDownPos.y - convey2DeviceDist, -1);
					Vector2 outLeftPos(outDeviceLeftPos.x - convey2DeviceDist, outDeviceLeftPos.y, -1);
					Vector2 outRightPos(outDeviceRightPos.x + convey2DeviceDist, outDeviceRightPos.y, -1);

					Vector2 inUpPos(inDeviceUpPos.x, inDeviceUpPos.y + convey2DeviceDist, -1);
					Vector2 inDownPos(inDeviceDownPos.x, inDeviceDownPos.y - convey2DeviceDist, -1);
					Vector2 inLeftPos(inDeviceLeftPos.x - convey2DeviceDist, inDeviceLeftPos.x, -1);
					Vector2 inRightPos(outDeviceRightPos.x + convey2DeviceDist, outDeviceRightPos.y, -1);


					//Ҫ�õ��ĸ��ֱȽ�ֵ
					double inR2outXDist, out2inLXDist, in2outLXDist, outR2inXDist;
					double outU2inDYDist, inR2outLXDist, outR2inLXDist;
					double outL2inLXDist, inU2outUYDist, inU2outDYDist, outU2inYDist;
					double out2inDYDist, inU2outYDist, in2outDYDist;
					//���ж����ߵķ�λ
					//��һ���ж϶��߷�λ�Ĳ���
					Vector2 inDevicePos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					Vector2 outDevicePos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					bool inLeft2OutBool, inRight2OutBool, inUp2OutBool, inDown2OutBool;

					int OutPosIndex_X, OutPosIndex_Y, InPosIndex_X, InPosIndex_Y;
					OutPosIndex_X = outDeviceIndex * 3; OutPosIndex_Y = outDeviceIndex * 3 + 1;
					InPosIndex_X = inDeviceIndex * 3; InPosIndex_Y = inDeviceIndex * 3 + 1;
					double moveLength = 0.0;
					switch (pointDirectArray[outPointDirect * 5 + inPointDirect])
					{
					case 1://����
					{
						//1.������y�Ϲ��ڽӽ�,��������y�϶���
						if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)
						{
							moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
							position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
							position_GPU[index * dim + InPosIndex_Y] += moveLength;
							break;
						}
						//2.����豸�ڳ����豸��/���Ͻ�
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						if (inUp2OutBool)
						{
							inR2outXDist = inRightPos.x - outPointTPos.x;
							if (inR2outXDist > 0 && inR2outXDist < conveyMinDist)
							{
								moveLength = inR2outXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							out2inLXDist = outPointTPos.x - inLeftPos.x;
							if (out2inLXDist > 0 && out2inLXDist < conveyMinDist)
							{
								moveLength = out2inLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						//3.����豸�ڳ����豸��/���½�
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							in2outLXDist = inPointTPos.x - outLeftPos.x;
							if (in2outLXDist > 0 && in2outLXDist < conveyMinDist)
							{
								moveLength = in2outLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							outR2inXDist = outRightPos.x - inPointTPos.x;
							if (outR2inXDist > 0 && outR2inXDist < conveyMinDist)
							{
								moveLength = outR2inXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						break;
					}
					case 11://����
					{
						//1.������y�Ϲ��ڽӽ�,��������y�϶���
						if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)//�������
						{
							moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
							position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
							position_GPU[index * dim + InPosIndex_Y] += moveLength;
							break;
						}
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						//2.����豸�ڳ����豸��/���Ͻ�
						if (inUp2OutBool)
						{
							in2outLXDist = inPointTPos.x - outLeftPos.x;
							if (in2outLXDist > 0 && in2outLXDist < conveyMinDist)
							{
								moveLength = in2outLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							outR2inXDist = outRightPos.x - inPointTPos.x;
							if (outR2inXDist > 0 && outR2inXDist < conveyMinDist)
							{
								moveLength = outR2inXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						//3.����豸�ڳ����豸��/���½�
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							inR2outXDist = inRightPos.x - outPointTPos.x;
							if (inR2outXDist > 0 && inR2outXDist < conveyMinDist)
							{
								moveLength = inR2outXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							out2inLXDist = outPointTPos.x - inLeftPos.x;
							if (out2inLXDist > 0 && out2inLXDist < conveyMinDist)
							{
								moveLength = out2inLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						break;
					}
					case 16://����
					{
						//1.������x�Ϲ��ڽӽ�,��������x�϶���
						if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
						{
							moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
							position_GPU[index * dim + OutPosIndex_X] -= moveLength;
							position_GPU[index * dim + InPosIndex_X] += moveLength;
							break;
						}
						//2.in��out���
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//in��out����
							out2inDYDist = outPointTPos.y - inDownPos.y;
							if (out2inDYDist > 0 && out2inDYDist < conveyMinDist)
							{
								moveLength = out2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in��out����
							inU2outYDist = inUpPos.y - outPointTPos.y;
							if (inU2outYDist > 0 && inU2outYDist < conveyMinDist)
							{
								moveLength = inU2outYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						//3.in��out�ұ�
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in��out����
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in��out����
							in2outDYDist = inPointTPos.y - outDownPos.y;
							if (in2outDYDist > 0 && in2outDYDist < conveyMinDist)
							{
								moveLength = in2outDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						break;
					}
					case 6://����
					{
						//1.������x�Ϲ��ڽӽ�,��������x�϶���
						if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
						{
							moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
							position_GPU[index * dim + OutPosIndex_X] -= moveLength;
							position_GPU[index * dim + InPosIndex_X] += moveLength;
							break;
						}
						//2.in��out�ұ�
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in��out����
							out2inDYDist = outPointTPos.y - inDownPos.y;
							if (out2inDYDist > 0 && out2inDYDist < conveyMinDist)
							{
								moveLength = out2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in��out����
							inU2outYDist = inUpPos.y - outPointTPos.y;
							if (inU2outYDist > 0 && inU2outYDist < conveyMinDist)
							{
								moveLength = inU2outYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						//3.in��out���
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//in��out����
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in��out����
							in2outDYDist = inPointTPos.y - outDownPos.y;
							if (in2outDYDist > 0 && in2outDYDist < conveyMinDist)
							{
								moveLength = in2outDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						break;
					}
					case 3://����
					{
						//1.in��out�Ϸ�������x�Ϻܽӽ�
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						if (inUp2OutBool)
						{
							if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
							{
								moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						//2.����豸�ڳ����豸��/���Ͻ�
						inUp2OutBool = IsInUp2Out(inDeviceDownPos, outDeviceUpPos);
						if (inUp2OutBool)//���ҿ�����һ�ַ�ʽ����
						{
							outU2inDYDist = outUpPos.y - inDownPos.y;
							if (outU2inDYDist > 0 && outU2inDYDist < conveyMinDist)
							{
								moveLength = outU2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//����豸�ڳ����豸����
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							//3.����豸�ڳ����豸����/��Ϊ����
							inR2outLXDist = inRightPos.x - outLeftPos.x;
							if (inR2outLXDist > 0 && inR2outLXDist < conveyMinDist)
							{
								moveLength = inR2outLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							outR2inLXDist = outRightPos.x - inLeftPos.x;
							if (outR2inLXDist > 0 && outR2inLXDist < conveyMinDist)
							{
								moveLength = outR2inLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						break;
					}
					case 14://����
					{
						//1.in��out��ߣ�����y�Ϻܽӽ�
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)
							{
								moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//2.in��out���
						inLeft2OutBool = IsInLeft2Out(inDeviceRightPos, outDeviceLeftPos);
						if (inLeft2OutBool)
						{
							inR2outLXDist = inRightPos.x - outLeftPos.x;
							if (inR2outLXDist > 0 && inR2outLXDist < conveyMinDist)
							{
								moveLength = inR2outLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
						}
						//3.in��out�ұߣ���Ϊ���Ϻ�����
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in��out����
							inU2outDYDist = inUpPos.y - outDownPos.y;
							if (inU2outDYDist > 0 && inU2outDYDist < conveyMinDist)
							{
								moveLength = inU2outDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							//in��out����
							outU2inDYDist = outUpPos.y - inDownPos.y;
							if (outU2inDYDist > 0 && outU2inDYDist < conveyMinDist)
							{
								moveLength = outU2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						break;
					}
					case 4://����(�Ȳ����ǿ��ܵ����)
					{
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//1.in��out����
							outU2inDYDist = outUpPos.y - inDownPos.y;
							if (outU2inDYDist > 0 && outU2inDYDist < conveyMinDist)
							{
								moveLength = outU2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						if (inUp2OutBool)
						{
							//2.in��out����
							out2inLXDist = outPointTPos.x - inLeftPos.x;
							if (out2inLXDist > 0 && out2inLXDist < conveyMinDist)
							{
								moveLength = out2inLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//3.in��out����һ��
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							//4.in��out���½�
							outR2inLXDist = outRightPos.x - inLeftPos.x;
							if (outR2inLXDist > 0 && outR2inLXDist < conveyMinDist)
							{
								moveLength = outR2inLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] -= moveLength;
								position_GPU[index * dim + InPosIndex_X] += moveLength;
								break;
							}
						}
						break;
					}
					case 2://����
					{
						//in��out��
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						if (inUp2OutBool)
						{
							//1.in��out����
							inR2outXDist = inRightPos.x - outPointTPos.x;
							if (inR2outXDist > 0 && inR2outXDist < conveyMinDist)
							{
								moveLength = inR2outXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
						}
						//in��out��
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);//�������Ҳ�����Ż�����ǰ���
						if (inRight2OutBool)
						{
							//2.in��out����
							outU2inDYDist = outUpPos.y - inDownPos.y;
							if (outU2inDYDist > 0 && outU2inDYDist < conveyMinDist)
							{
								moveLength = outU2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//in��out��
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//3.in��out����
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//in��out��
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							//3.in��out����
							inR2outLXDist = inRightPos.x - outLeftPos.x;
							if (inR2outLXDist > 0 && inR2outLXDist < conveyMinDist)
							{
								moveLength = inR2outLXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
						}
						break;
					}
					default:
						break;
					}
				}
				curDeviceLink_Index++;
			}
		}
#pragma endregion

#pragma region �������ڵ�ļ�������
		//����in��out�㱻�ֿ���
		//Ȼ������е�inoutPoint��ֵ
		//InoutPoint* tempInoutPoints = new InoutPoint[totalInPoint + totalOutPoint];
		int ioPIndex = 0;
		for (int i = 0; i < DeviceSum; i++)
		{
			for (int pointIndex = 0; pointIndex < adjPInCount[i]; ++pointIndex)
			{
				AdjPoint& point = adjPointsIn[accumAdjPInCount[i] + pointIndex];
				InoutPoint ioPoint;//
				ioPoint.pointDirect = point.direct;
				Vector2 axis(point.pos.x + position_GPU[index * dim + 3 * i], point.pos.y + position_GPU[index * dim + 3 * i + 1], -1);
				ioPoint.pointAxis = axis;
				inoutPoints[index * inoutPSize + ioPIndex] = ioPoint;//ע��ƫ��ֵ
				ioPIndex++;
			}
			for (int pointIndex = 0; pointIndex < adjPOutCount[i]; ++pointIndex)
			{
				AdjPoint& point = adjPointsOut[accumAdjPOutCount[i] + pointIndex];
				InoutPoint ioPoint;
				ioPoint.pointDirect = point.direct;
				Vector2 axis(point.pos.x + position_GPU[index * dim + 3 * i], point.pos.y + position_GPU[index * dim + 3 * i + 1], -1);
				ioPoint.pointAxis = axis;
				inoutPoints[index * inoutPSize + ioPIndex] = ioPoint;//ע��ƫ��ֵ
				ioPIndex++;
			}
		}

#pragma endregion

#pragma region �����豸����ͳ�������깹��·����ͼ
		double* horizonAxisList = new double[horiPointCount];
		double* verticalAxisList = new double[vertPointCount];
		int curHoriIndex = 0;
		int curVertIndex = 0;
		//�ȶԳ���ڵ�ˮƽ�ʹ�ֱ���з���(ע�����ƫ����)
		//�����ˮƽ�ʹ�ֱ����ڵ����Ŀ��һ����horiCount��vertCount��
		for (int i = 0; i < DeviceSum; i++)
		{
			for (int pointIndex = 0; pointIndex < adjPInCount[i]; pointIndex++)
			{
				AdjPoint& p = adjPointsIn[accumAdjPInCount[i] + pointIndex];
				if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//����
				{
					horizonAxisList[curHoriIndex++] = p.pos.x + position_GPU[index * dim + i * 3];
				}
				else {//����
					verticalAxisList[curVertIndex++] = p.pos.y + position_GPU[index * dim + i * 3 + 1];
				}
			}
			for (int pointIndex = 0; pointIndex < adjPOutCount[i]; pointIndex++)
			{
				AdjPoint& p = adjPointsOut[accumAdjPOutCount[i] + pointIndex];
				if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//����
				{
					horizonAxisList[curHoriIndex++] = p.pos.x + position_GPU[index * dim + i * 3];
				}
				else {//����
					verticalAxisList[curVertIndex++] = p.pos.y + position_GPU[index * dim + i * 3 + 1];
				}
			}
		}
		//�ֿ����2����&���ڣ�horiCount��vertCount+=2��
		horizonAxisList[curHoriIndex++] = entrancePos.x;
		verticalAxisList[curVertIndex++] = entrancePos.y;

		horizonAxisList[curHoriIndex++] = exitPos.x;
		verticalAxisList[curVertIndex++] = exitPos.y;
		//����ÿ���豸������ĸ���Χ����Ϊ�����ϰ���ķ�Χ��
		double* DeviceLowXList = new double[DeviceSum];
		double* DeviceHighXList = new double[DeviceSum];
		double* DeviceLowYList = new double[DeviceSum];
		double* DeviceHighYList = new double[DeviceSum];
		//ÿ���豸��Χ��4����
		for (int i = 0; i < dim; i += 3) {
			outSizeLength = 0.5 * size_copy[i / 3].x + convey2DeviceDist;
			outSizeWidth = 0.5 * size_copy[i / 3].y + convey2DeviceDist;
			double LowX = position_GPU[index * dim + i] - outSizeLength;
			double HighX = position_GPU[index * dim + i] + outSizeLength;
			double LowY = position_GPU[index * dim + i + 1] - outSizeWidth;
			double HighY = position_GPU[index * dim + i + 1] + outSizeWidth;

			verticalAxisList[curVertIndex++] = LowY;
			verticalAxisList[curVertIndex++] = HighY;
			horizonAxisList[curHoriIndex++] = LowX;
			horizonAxisList[curHoriIndex++] = HighX;

			//ÿ���豸���ĸ���Χ
			DeviceLowXList[i / 3] = LowX;
			DeviceHighXList[i / 3] = HighX;
			DeviceLowYList[i / 3] = LowY;
			DeviceHighYList[i / 3] = HighY;

			//��ֹ·�������豸�ڲ�
			horizonAxisList[curHoriIndex++] = position_GPU[index * dim + i];
			verticalAxisList[curVertIndex++] = position_GPU[index * dim + i + 1];

		}
		//��һ��������Щ��
		//����Щ������갴�մ�С��������
		//sort������Ҫ�Լ�ʵ�֣�����ʹ��cuda�⺯����
		Double_Sort(horizonAxisList, 0, horiPointCount - 1);
		Double_Sort(verticalAxisList, 0, vertPointCount - 1);
		//ֻ�������ظ��ĵ㣨�Լ�ʵ�֣�

		int uniqueHoriPCount = horiPointCount;
		int uniqueVertPCount = vertPointCount;
		int unique_end1 = Double_Unique(horizonAxisList, 0, uniqueHoriPCount - 1);
		uniqueHoriPCount = unique_end1;
		//verticalAxisList.erase(unique_end1, verticalAxisList.end());
		int unique_end2 = Double_Unique(verticalAxisList, 0, uniqueVertPCount - 1);
		uniqueVertPCount = unique_end2;
		//horizonAxisList.erase(unique_end2, horizonAxisList.end());

		//�����е��ϰ�����±�
		int barrierRowNum = 200;//��һ���̶���С�����ڴ�
		int* barrierRowIndexList = new int[barrierRowNum];
		int barrierColNum = 200;
		int* barrierColIndexList = new int[barrierColNum];
		int totalBarRowNum = 0;//��¼barrier��ʵ������Ŀ
		int totalBarColNum = 0;//��¼barrier��ʵ������Ŀ

		//����Щ����ȥ���·����map��map�Ƕ�ά�ģ��൱�ڶ�ά�����
		//��һά�����ά
		//horiNum��ӦcolNum��vertNum��ӦrowNum
		int pathColNum = uniqueHoriPCount;
		int pathRowNum = uniqueVertPCount;
		APoint** pathPointMap = new APoint*[pathColNum * pathRowNum];

		for (int i = 0; i < pathColNum * pathRowNum; i++)
		{
			pathPointMap[i] = new APoint();//���������ڴ�
			int rowIndex = i / pathColNum;//��Ӧ�����i
			int colIndex = i % pathColNum;//��Ӧ�����j
			pathPointMap[i]->x = horizonAxisList[colIndex];
			pathPointMap[i]->y = verticalAxisList[rowIndex];
			//���������豸�����Ƿ��к�������ص��ģ�ʵ�ֱ���ϰ��㣩
			for (int k = 0; k < DeviceSum; k++)
			{
				if (pathPointMap[i]->x - DeviceLowXList[k] >= 0.01 && DeviceHighXList[k] - pathPointMap[i]->x >= 0.01
					&& pathPointMap[i]->y - DeviceLowYList[k] >= 0.01 && DeviceHighYList[k] - pathPointMap[i]->y >= 0.01)
				{
					pathPointMap[i]->type = AType::ATYPE_BARRIER;
					//�ϰ�����ͼ�е��±�
					barrierRowIndexList[totalBarRowNum++] = rowIndex;
					barrierColIndexList[totalBarColNum++] = colIndex;
				}
			}
			pathPointMap[i]->colIndex = colIndex;//��¼�µ���·����map�е��±�
			pathPointMap[i]->rowIndex = rowIndex;
		}
#pragma endregion

#pragma region Ѱ·
		CAstar* star = new CAstar();
		//ֻ��kernal�����ڲ�ִ�У���new��delete����
		star->_allPoints = pathPointMap;
		star->pointColNum = pathColNum;
		star->pointRowNum = pathRowNum;

		int beginRowIndex, beginColIndex, endRowIndex, endColIndex;

		double totalTime = 0.0;
		//�����е�·����Ϣ
		PointLink* copyPLinks = new PointLink[totalLinkSum];//��ʱÿ��·����50����
		int curLinkIndex = 0;//��ǰ��link�±�
		int totalLinkPointSum = 0;//����link�����е����Ŀ��ÿ��link�кܶ�㣩
		for (int i = 0; i < CargoTypeNum; i++)//��������
		{
			//CargoType curCargoType = proParas.cargoTypeList[i];
			//ÿ�ֻ�����ܾ�������豸��Ҳ����һ��Type��Ӧ���link
			for (int j = 0; j < linkSum[i]; j++)
			{
				PathDirection pathBeginDirect;
				PathDirection pathEndDirect;
				int forwardDeviceIndex, curDeviceIndex;//�豸1���豸2ID
				int forwardOutIndex, curInIndex;//����ڵ��±�
				double device1PosX, device1PosY, device2PosX, device2PosY;//�豸��Χ���ĸ���
				double initDevice1PosX, initDevice1PosY, initDevice2PosX, initDevice2PosY;//����δ���Ӱ�Χ�ߵ�����

				double deviceDistance = 0.0;//����
				double outDSizeL, inDSizeL;

				forwardDeviceIndex = deviceLinkList[accumLinkSum[i] + j].outDeviceIndex;
				curDeviceIndex = deviceLinkList[accumLinkSum[i] + j].inDeviceIndex;

				forwardOutIndex = deviceLinkList[accumLinkSum[i] + j].outPointIndex;
				curInIndex = deviceLinkList[accumLinkSum[i] + j].inPointIndex;
				if (forwardDeviceIndex == -1)//˵���ǲֿ����
				{
					/*forwardOutIndex = 0;
					curInIndex = proParas.cargoTypeList[i].deviceLinkList[j].inPointIndex;*/

					//��ͷ�ͽ�β��ĳ�����
					pathBeginDirect = PathDirection::Vertical;
					pathEndDirect = (adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].direct % 2 == 0)
						? PathDirection::Horizon : PathDirection::Vertical;
					outDSizeL = 0.0;
					inDSizeL = pathEndDirect == PathDirection::Horizon ? (0.5 * size_copy[curDeviceIndex].x) : (0.5 * size_copy[curDeviceIndex].y);

					device1PosX = entrancePos.x;
					device1PosY = entrancePos.y;
					device2PosX = adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].pos.x + position_GPU[index * dim + curDeviceIndex * 3];
					device2PosY = adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].pos.y + position_GPU[index * dim + curDeviceIndex * 3 + 1];

					initDevice1PosX = device1PosX;
					initDevice1PosY = device1PosY;
					initDevice2PosX = device2PosX;
					initDevice2PosY = device2PosY;
					//�õ��豸��Χ�ĵ�
					switch (adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].direct)
					{
					case PointDirect::Up:
						device2PosY += convey2DeviceDist;
						break;
					case PointDirect::Down:
						device2PosY -= convey2DeviceDist;
						break;
					case PointDirect::Left:
						device2PosX -= convey2DeviceDist;
						break;
					case PointDirect::Right:
						device2PosX += convey2DeviceDist;
						break;
					}
				}
				else if (curDeviceIndex == -2)//˵���ǳ���
				{
					//��ͷ�ͽ�β��ĳ����� 
					pathBeginDirect = (adjPointsOut[accumAdjPInCount[forwardDeviceIndex] + forwardOutIndex].direct % 2 == 0)
						? PathDirection::Horizon : PathDirection::Vertical;
					pathEndDirect = PathDirection::Vertical;

					outDSizeL = pathBeginDirect == PathDirection::Horizon ? (0.5 * size_copy[forwardDeviceIndex].x) : (0.5 * size_copy[forwardDeviceIndex].y);
					inDSizeL = 0.0;

					device1PosX = adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].pos.x + position_GPU[index * dim + forwardDeviceIndex * 3];
					device1PosY = adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].pos.y + position_GPU[index * dim + forwardDeviceIndex * 3 + 1];

					device2PosX = exitPos.x;
					device2PosY = exitPos.y;

					initDevice1PosX = device1PosX;
					initDevice1PosY = device1PosY;
					initDevice2PosX = device2PosX;
					initDevice2PosY = device2PosY;
					//�õ��豸��Χ�ĵ� 
					switch (adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].direct)
					{
					case PointDirect::Up:
						//device1PosY += copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosY += convey2DeviceDist;
						break;
					case PointDirect::Down:
						//device1PosY -= copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosY -= convey2DeviceDist;
						break;
					case PointDirect::Left:
						//device1PosX -= copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosX -= convey2DeviceDist;
						break;
					case PointDirect::Right:
						//device1PosX += copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosX += convey2DeviceDist;
						break;
					}
				}
				else//��ͨ
				{
					//forwardOutIndex = proParas.cargoTypeList[i].deviceLinkList[j].outPointIndex - 1;
					//curInIndex = proParas.cargoTypeList[i].deviceLinkList[j].inPointIndex - 1;

					//��ͷ�ͽ�β��ĳ�����
					pathBeginDirect = (adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].direct % 2 == 0)
						? PathDirection::Horizon : PathDirection::Vertical;
					pathEndDirect = (adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].direct % 2 == 0)
						? PathDirection::Horizon : PathDirection::Vertical;


					outDSizeL = pathBeginDirect == PathDirection::Horizon ? (0.5 * size_copy[forwardDeviceIndex].x) : (0.5 * size_copy[forwardDeviceIndex].y);
					inDSizeL = pathEndDirect == PathDirection::Horizon ? (0.5 * size_copy[curDeviceIndex].x) : (0.5 * size_copy[curDeviceIndex].y);

					device1PosX = adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].pos.x + position_GPU[index * dim + forwardDeviceIndex * 3];
					device1PosY = adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].pos.y + position_GPU[index * dim + forwardDeviceIndex * 3 + 1];
					device2PosX = adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].pos.x + position_GPU[index * dim + curDeviceIndex * 3];
					device2PosY = adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].pos.y + position_GPU[index * dim + curDeviceIndex * 3 + 1];


					initDevice1PosX = device1PosX;
					initDevice1PosY = device1PosY;
					initDevice2PosX = device2PosX;
					initDevice2PosY = device2PosY;
					//�õ��豸��Χ�ĵ� 
					switch (adjPointsOut[accumAdjPOutCount[forwardDeviceIndex] + forwardOutIndex].direct)
					{
					case PointDirect::Up:
						//device1PosY += copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosY += convey2DeviceDist;
						break;
					case PointDirect::Down:
						//device1PosY -= copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosY -= convey2DeviceDist;
						break;
					case PointDirect::Left:
						//device1PosX -= copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosX -= convey2DeviceDist;
						break;
					case PointDirect::Right:
						//device1PosX += copyDeviceParas[forwardDeviceIndex].spaceLength;
						device1PosX += convey2DeviceDist;
						break;
					}
					switch (adjPointsIn[accumAdjPInCount[curDeviceIndex] + curInIndex].direct)
					{
					case PointDirect::Up:
						//device2PosY += copyDeviceParas[curDeviceIndex].spaceLength;
						device2PosY += convey2DeviceDist;
						break;
					case PointDirect::Down:
						//device2PosY -= copyDeviceParas[curDeviceIndex].spaceLength;
						device2PosY -= convey2DeviceDist;
						break;
					case PointDirect::Left:
						//device2PosX -= copyDeviceParas[curDeviceIndex].spaceLength;
						device2PosX -= convey2DeviceDist;
						break;
					case PointDirect::Right:
						//device2PosX += copyDeviceParas[curDeviceIndex].spaceLength;
						device2PosX += convey2DeviceDist;
						break;
					}
				}
				//�������·��
				beginRowIndex = FindAxisIndex(device1PosY, verticalAxisList, uniqueVertPCount);
				beginColIndex = FindAxisIndex(device1PosX, horizonAxisList, uniqueHoriPCount);
				endRowIndex = FindAxisIndex(device2PosY, verticalAxisList, uniqueVertPCount);
				endColIndex = FindAxisIndex(device2PosX, horizonAxisList, uniqueHoriPCount);

				//�õ�·����path�ǵ�һ���ڵ�
				APoint* path = findWay(star->curPathDirect, pathBeginDirect, star->_allPoints, star->_endPoint, star->_neighbourList,
					star->_curPoint, star->_openList, star->_closeList, star->openList_CurSize, star->closeList_CurSize, star->neighbourList_CurSize,
					star->pointColNum, star->pointRowNum, beginRowIndex, beginColIndex, endRowIndex, endColIndex);
				//�����еĽ⣬ֱ���˳�
				if (path == nullptr)
				{
					fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
					printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
					return;
				}
				//����·�����㳤��
				deviceDistance = CalcuPathLength(path) + outDSizeL + inDSizeL;


#pragma region ����·����ֻ������ʼ�� + ·���е�ת��㣩
				//·����������
				//���·���еĵ㻹û�м򻯣��ȸ�ÿ��·�����ô�СΪ50
				Vector2* points1 = new Vector2[fixedLinkPointSum];
				int points1Index = 0;
				Vector2 endP1(initDevice2PosX, initDevice2PosY, -1);
				points1[points1Index++] = endP1;
				APoint* copyPath = path;
				while (copyPath)
				{
					Vector2 tempP(copyPath->x, copyPath->y, -1);
					points1[points1Index++] = tempP;
					copyPath = copyPath->parent;
				}
				Vector2 startP1(initDevice1PosX, initDevice1PosY, -1);
				points1[points1Index++] = startP1;//points1Index��ĳһ��link��ʵ�ʵ����Ŀ
				totalLinkPointSum += points1Index;//ͳ��link�����е�ʵ�ʵ����Ŀ������֮�����ռ䣩

				PointLink pointLink1(forwardDeviceIndex, forwardOutIndex, curDeviceIndex, curInIndex, points1, points1Index, -1);
				copyPLinks[curLinkIndex++] = pointLink1;

				//Vector2 lastP(path->x, path->y);//�����һ���㿪ʼ
				//points.push_back(lastP);
				//path = path->parent;
				//while (path)
				//{
				//	Vector2 curP(path->x, path->y);
				//	//ֻ��ת��ĵ�Żᱻ����·��
				//	if (pathCurDirect == PathDirection::Horizon && curP.x == lastP.x)
				//	{
				//		//if (curIterNum != 0)
				//			//if ((lastP.x != points.back().x || lastP.y != points.back().y)
				//			//	&& CalcuDeviceDist(lastP, points.back()) < proParas.conveyMinDist)
				//			//{
				//			//	//punishValue1 = 200 * (curIterNum + 1);
				//			//	//punishValue2 = 40 * (curIterNum + 1);
				//			//	particle.fitness_[0] = particle.fitness_[1] = MAX_FITNESS;
				//			//	return;
				//			//}
				//		points.push_back(lastP);
				//		pathCurDirect = PathDirection::Vertical;
				//	}
				//	else if (pathCurDirect == PathDirection::Vertical && curP.y == lastP.y)
				//	{
				//		//if (curIterNum != 0) {
				//			//if ((lastP.x != points.back().x || lastP.y != points.back().y)
				//			//	&& CalcuDeviceDist(lastP, points.back()) < proParas.conveyMinDist)
				//			//{
				//			//	//punishValue1 = 200 * (curIterNum + 1);
				//			//	//punishValue2 = 40 * (curIterNum + 1);
				//			//	particle.fitness_[0] = particle.fitness_[1] = MAX_FITNESS;
				//			//	return;
				//			//}
				//		points.push_back(lastP);
				//		pathCurDirect = PathDirection::Horizon;
				//	}
				//	path = path->parent;
				//	lastP = curP;
				//}
				////if (curIterNum != 0)
				//	//if (CalcuDeviceDist(lastP, points.back()) < proParas.conveyMinDist)
				//	//{
				//	//	//punishValue1 = 200 * (curIterNum + 1);
				//	//	//punishValue2 = 40 * (curIterNum + 1);
				//	//	particle.fitness_[0] = particle.fitness_[1] = MAX_FITNESS;
				//	//	return;
				//	//}
				//points.push_back(lastP);
				//Vector2 startP(initDevice1PosX, initDevice1PosY);
				//points.push_back(startP);

				//PointLink pointLink(forwardDeviceIndex, forwardOutIndex, curDeviceIndex, curInIndex, points);
				//particle.pointLinks.push_back(pointLink);
#pragma endregion

				//��������ʱ��(�������� * ·�߳��� * ����Ч��) 
				totalTime += totalVolume[i] * deviceDistance * conveySpeed;
				//�����豸����ʱ��(�������� * ����Ч��)
				//totalTime += curCargoType.totalVolume * curDevice.workSpeed;

				resetAStar(star->pointRowNum, star->pointColNum, star->_allPoints, star->_openList, star->_closeList, star->_neighbourList,
					star->openList_CurSize, star->closeList_CurSize, star->neighbourList_CurSize);
				//���ϰ������±��
				for (int i = 0; i < totalBarRowNum; i++)
				{
					star->_allPoints[barrierRowIndexList[i] * pathColNum + barrierColIndexList[i]]->type = AType::ATYPE_BARRIER;
				}
			}
		}
#pragma endregion

#pragma region �����ֽ��ת��Ϊ���ͻ�����(ͬʱ��ǿһ����������̾���Լ��)

		//ֻ����ת��������
		int segPathSet_PointSum = fixedLinkPointSum * totalLinkSum;//50����
		int segPathSet_CurIndex = 0;
		SegPath* segPathSet = new SegPath[segPathSet_PointSum];

		//////���������ظ���·�߶Σ���ȡsegPathSet//////
		//1.��copyPLinks������е�ŵ�seg������
		for (int i = 0; i < totalLinkSum; i++)
		{
			for (int j = copyPLinks[i].pointNum - 1; j > 0; j--)
			{
				Vector2Int p1 = Multi10000ToInt(copyPLinks[i].points[j]);
				Vector2Int p2 = Multi10000ToInt(copyPLinks[i].points[j - 1]);
				if (p1.ANEqualB(p2, -1))//���겻��һ��
				{
					SegPath temp(p1, p2, -1);
					segPathSet[segPathSet_CurIndex++] = temp;
				}
			}
		}
		//2.��seg�����������+ȥ��
		SegPath_Sort(segPathSet, 0, segPathSet_CurIndex - 1);
		//ȥ�غ�ĵ����Ŀ
		int segPathSet_UniqueSum = SegPath_Unique(segPathSet, 0, segPathSet_CurIndex - 1);

		//////����set�������set�����е� ������·�ߵĴ�ֱˮƽ��Ŀ�����Ƴ���ȣ�&�Ƿ��� ��Ϣ//////
		//��Щ��Ϣ�浽pathPointInfoMap��
		//map<Vector2Int, PointInfo> pathPointInfoMap;


		//·������Ϣmap(ÿ���������:ÿ�������Ϣ)
		int pathPointInfoMap_PointSum = segPathSet_UniqueSum * 2;//ÿ��seg��Ӧ������
		PointInfo* pathPointInfoMap = new PointInfo[pathPointInfoMap_PointSum];
		int pathPointInfoMap_CurIndex = 0;//��ǰ�±�

		//ͬ���ģ�ֻ�����������map
		//���ʵ�֣�
		//1.�ֽ�segPath�����еĵ�ŵ�PointInfo��
		for (int segIndex = 0; segIndex < segPathSet_UniqueSum; segIndex++)
		{
			if (segPathSet[segIndex].direct == PathPointDirect::Vert) {
				pathPointInfoMap[pathPointInfoMap_CurIndex++] = PointInfo(segPathSet[segIndex].p1, 1, 0, false, -1);
				pathPointInfoMap[pathPointInfoMap_CurIndex++] = PointInfo(segPathSet[segIndex].p2, 1, 0, false, -1);
			}
			else
			{
				pathPointInfoMap[pathPointInfoMap_CurIndex++] = PointInfo(segPathSet[segIndex].p1, 0, 1, false, -1);
				pathPointInfoMap[pathPointInfoMap_CurIndex++] = PointInfo(segPathSet[segIndex].p2, 0, 1, false, -1);
			}
		}
		//2.�����е��������Ŀ���ǽ�������ͬ�ĵ�ۼ���һ��
		PointInfo_Sort(pathPointInfoMap, 0, pathPointInfoMap_CurIndex - 1);
		//3.������������飬����ÿ�����vertNum��horiNum��Ŀ,˳��ȥ��
		int pathPointInfoMap_UniqueSum = PointInfo_CalcuAndUnique(pathPointInfoMap, 0, pathPointInfoMap_CurIndex - 1);


		//4.���ݵ��vertNum��horiNum���ж�ÿ�����Ƿ񱻱���
		for (int i = 0; i < pathPointInfoMap_UniqueSum; i++)
		{
			if ((pathPointInfoMap[i].horiDirNum == 1 && pathPointInfoMap[i].vertDirNum == 0)
				|| (pathPointInfoMap[i].horiDirNum == 0 && pathPointInfoMap[i].vertDirNum == 1)) {
				pathPointInfoMap[i].isKeep = true;
			}
			if (pathPointInfoMap[i].horiDirNum >= 1 && pathPointInfoMap[i].vertDirNum >= 1) {
				pathPointInfoMap[i].isKeep = true;
			}
		}

		//��pathPointInfoMap�еõ�tempStrConveyorList��tempCurveConveyorList
		//ԭ��������set��map������ʵ�ֿ������
		//����û���ˣ�ֻ���ö��ֲ�����

		//ֱ�����ͻ���Ϣ�б�
		int tempStrConveyorList_PointSum = fixedUniqueLinkPointSum * totalLinkSum;//20����
		//StraightConveyorInfo* tempStrConveyorList = new StraightConveyorInfo[tempStrConveyorList_PointSum];
		int tempStrConveyorList_CurIndex = 0;

		//ת�����ͻ���Ϣ�б�
		int tempCurveConveyorList_PointSum = fixedUniqueLinkPointSum * totalLinkSum;//20����
		//Vector2Int* tempCurveConveyorList = new Vector2Int[tempCurveConveyorList_PointSum];
		int tempCurveConveyorList_CurIndex = 0;

		for (int i = 0; i < totalLinkSum; i++) {
			StraightConveyorInfo tempStrInfo;
			tempStrInfo.startPos = Multi10000ToInt(copyPLinks[i].points[copyPLinks[i].pointNum - 1]);//��ͷ
			PointInfo curPointInfo = FindPointInfo(pathPointInfoMap, 0, pathPointInfoMap_UniqueSum - 1, tempStrInfo.startPos);
			tempStrInfo.startVnum = curPointInfo.vertDirNum;
			tempStrInfo.startHnum = curPointInfo.horiDirNum;
			for (int j = copyPLinks[i].pointNum - 2; j >= 0; j--) {
				Vector2Int p = Multi10000ToInt(copyPLinks[i].points[j]);
				PointInfo pPathPoint = FindPointInfo(pathPointInfoMap, 0, pathPointInfoMap_UniqueSum - 1, p);
				if (pPathPoint.isKeep == true) {//���������ظ���
					//�ȸ���ֱ�����ͻ�
					if (tempStrInfo.startPos.ANEqualB(p, -1)) {
						tempStrInfo.endPos = p;
						//if (curIterNum != 0)
						if (tempStrInfo.startPos.Distance(tempStrInfo.endPos) < conveyMinDist)
						{
							//cout << "������̫��" << endl;
							//punishValue1 = 150 * (curIterNum + 1);
							//punishValue2 = 30 * (curIterNum + 1);
							fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
							printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
							return;
						}
						tempStrInfo.endVnum = pPathPoint.vertDirNum;
						tempStrInfo.endHnum = pPathPoint.horiDirNum;
						//����ط�����ȥ�ص�
						//���û��set��ֻ��������+ȥ�ظ���
						//ע��ÿ��p�����յ�Ҳ����һ�������
						strConveyorList[index * tempStrConveyorList_PointSum + tempStrConveyorList_CurIndex++] = tempStrInfo;//ע��ƫ��ֵ
						tempStrInfo.startPos = p;
						tempStrInfo.startVnum = pPathPoint.vertDirNum;
						tempStrInfo.startHnum = pPathPoint.horiDirNum;
					}
					//ֻҪ����ʼ�յ㣬����Ҫ����ת�����ͻ�
					if (!(pPathPoint.horiDirNum == 1 && pPathPoint.vertDirNum == 0)
						&& !(pPathPoint.horiDirNum == 0 && pPathPoint.vertDirNum == 1)) {
						//tempCurveConveyorList.insert(p);
						curveConveyorList[index * tempCurveConveyorList_PointSum + tempCurveConveyorList_CurIndex++] = p;
					}
				}
			}
		}
#pragma endregion

#pragma region ����Ŀ�꺯��ֵ
		//����ֱ�ߺ�ת�����ͻ���set���õ����ͻ����ܳɱ�
		double conveyorTotalCost = 0.0;
		for (int i = 0; i < tempStrConveyorList_CurIndex; i++) {
			conveyorTotalCost += strConveyorUnitCost * strConveyorList[index * tempStrConveyorList_PointSum + i]
				.startPos.Distance(strConveyorList[index * tempStrConveyorList_PointSum + i].endPos);
		}
		conveyorTotalCost += curveConveyorUnitCost * tempStrConveyorList_CurIndex;
		//������Ӧ��ֵ
		fitness_GPU[index * fitnessCount + 0] = totalTime;
		fitness_GPU[index * fitnessCount + 1] = conveyorTotalCost;

		//�����൱�ڼ�����˿��ܵ���������ߣ�������Ҫ����bestPathInfoList
		//������Ӧ���Ƿ�����ѡ�����BestPathInfoList
		//һ�����Բ��м������ֵ�Ĺ�Լ�㷨
		const int ParticleNum = 100;
		__shared__ int particleIndexList[ParticleNum];
		//Ϊ�丳ֵ
		particleIndexList[index] = index;
		//�ȴ�ÿ���̸߳�ֵ���
		__syncthreads();
		//ʵ�ֹ�Լ����������������������ӵ�index(���index�浽bestParticleIndex[0])
		int leng = particleNum;
		for (int i = particleNum / 2.0 + 0.5; i > 1; i = i / 2.0 + 0.5)
		{
			if (index < i)
			{
				if (index + i < leng)
				{
					//Ĭ�ϴ�����Ӧ��1
					if (fitness_GPU[particleIndexList[index] * fitnessCount + 0]
						< fitness_GPU[(particleIndexList[index] + i) * fitnessCount + 0]) {
						//����
						int temp = particleIndexList[index];
						particleIndexList[index] = particleIndexList[index + i];
						particleIndexList[index + i] = temp;
					}
				}
			}
			__syncthreads();
			leng = leng / 2.0 + 0.5;
		}
		if (index == 0)
		{
			bestParticleIndex[0] = fitness_GPU[particleIndexList[0] * fitnessCount + 0] > fitness_GPU[particleIndexList[1] * fitnessCount + 0]
				? particleIndexList[0] : particleIndexList[1];
		}
		__syncthreads();
		//bestParticleIndex[0]�������ѵ�fitness��Ӧ�������±�
		//for (int i = 0; i < fitnessCount; ++i) {
		//	if (fitness_GPU[index * fitnessCount + i] < bestPathInfoList[i].curBestFitnessVal) {//��Ҫ����
		//		bestPathInfoList[i].curBestFitnessVal = fitness_GPU[index * fitnessCount + i];
		//		bestPathInfoList[i].inoutPoints = inoutPoints + index * inoutPSize;//ע��ƫ��ֵ
		//		bestPathInfoList[i].inoutPSize = totalInPoint + totalOutPoint;//����������
		//		bestPathInfoList[i].strConveyorList = strConveyorList + index * tempStrConveyorList_PointSum;
		//		bestPathInfoList[i].strConveyorListSum = tempStrConveyorList_CurIndex;//��СҪ��¼����ʼ�����С��fixedUniqueLinkPointSum * totalLinkSum * fitnessCount
		//		bestPathInfoList[i].curveConveyorList = curveConveyorList + index * tempCurveConveyorList_PointSum;
		//		bestPathInfoList[i].curveConveyorListSum = tempCurveConveyorList_CurIndex;//fixedUniqueLinkPointSum * totalLinkSum * fitnessCount
		//	}
		//}
#pragma endregion
		delete[] DeviceLowXList;
		delete[] DeviceHighXList;
		delete[] DeviceLowYList;
		delete[] DeviceHighYList;

	}
#pragma endregion
	printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
	return;
}
//˳ʱ����ת�������
static __device__ Vector2 Rotate(Vector2 pointPos, Vector2 centerPos, float rotateAngle)
{
	float xx = (pointPos.x - centerPos.x) * cos(rotateAngle * (PI / 180)) + (pointPos.y - centerPos.y) * sin(rotateAngle * (PI / 180)) + centerPos.x;
	float yy = -(pointPos.x - centerPos.x) * sin(rotateAngle * (PI / 180)) + (pointPos.y - centerPos.y) * cos(rotateAngle * (PI / 180)) + centerPos.y;
	Vector2 result(xx, yy, -1);
	return result;
}
static __device__ int FindAxisIndex(double axis, const double* axisList, int axisCount)
{
	//�ö��ַ�����
	int low = 0;
	int high = axisCount - 1;
	int result = 0;
	while (low <= high)
	{
		int middle = (low + high) >> 1;
		if (abs(axisList[middle] - axis) <= 0.0001)
		{
			result = middle;
			break;
		}
		else if (axisList[middle] > axis)
		{
			high = middle - 1;
		}
		else if (axisList[middle] < axis)
		{
			low = middle + 1;
		}
	}
	return result;
}
////����ռ�����
//double CalcuTotalArea(Particle& particle, DevicePara* copyDeviceParas) {
//	double area = 0;
//	double min_X, min_Y, max_X, max_Y;
//	min_X = min_Y = INT_MAX;
//	max_X = max_Y = -INT_MAX;
//	for (int i = 0; i < particle.dim_; i += 3) {
//		double outSizeLength = copyDeviceParas[i / 3].size.x * 0.5 + copyDeviceParas[i / 3].spaceLength;
//		double outSizeWidth = copyDeviceParas[i / 3].size.y * 0.5 + copyDeviceParas[i / 3].spaceLength;
//		min_X = min(min_X, particle.position_[i] - outSizeLength);
//		max_X = max(max_X, particle.position_[i] + outSizeLength);
//		min_Y = min(min_Y, particle.position_[i + 1] - outSizeWidth);
//		max_Y = max(max_Y, particle.position_[i + 1] + outSizeWidth);
//	}
//	//���������
//	area = (max_X - min_X) * (max_Y - min_Y);
//	return area;
//}
//�ȳ���10000��Ȼ���������뵽Int
static __device__ int Multi10000ToInt(double num)
{
	//ʹ���Զ����Round����
	return MyRound(num * 10000);
}
static __device__ Vector2Int Multi10000ToInt(Vector2 v)
{
	return Vector2Int(Multi10000ToInt(v.x), Multi10000ToInt(v.y), -1);
}