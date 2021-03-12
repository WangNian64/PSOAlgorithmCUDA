#pragma once
#include <algorithm>
#include <math.h> 
#include "Tools.h"
#include "AStar.h"
#include "PSO.h"
#pragma region 判断两个坐标的上下或者左右关系
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
	/*ProblemParas proParas, 固定参数的，不用管*/
	int DeviceSum, int fixedLinkPointSum, int fixedUniqueLinkPointSum, int vertPointCount, int horiPointCount, double workShopLength, double workShopWidth, double convey2DeviceDist, /*double conveyWidth, */
	double strConveyorUnitCost, double curveConveyorUnitCost, double conveyMinDist, /*double conveyMinLength, */double conveySpeed, Vector2 entrancePos, Vector2 exitPos,
	int CargoTypeNum, int totalLinkSum,

	/*CargoType* 固定参数*/
	/*int* deviceSum, */int* linkSum, int* accumLinkSum, DeviceLink* deviceLinkList, double* totalVolume,

	/*DevicePara**/
	Vector2* size, double* spaceLength, int* adjPInCount, int* adjPOutCount, int* accumAdjPInCount, int* accumAdjPOutCount,
	int totalInPoint, int totalOutPoint, AdjPoint* adjPointsIn, AdjPoint* adjPointsOut,
	/*Particle*/
	int dim, int fitnessCount, double* fitness_GPU, double* position_GPU, /*double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU*/
	/*存储所有粒子输送线路信息*/
	double* curBestFitnessVal, int inoutPSize, InoutPoint* inoutPoints, StraightConveyorInfo* strConveyorList,
	int* strConveyorListSum, Vector2Int* curveConveyorList, int* curveConveyorListSum,
	int* pointDirectArray, curandState* globalState);

double CalcuTotalArea(Particle& particle, DevicePara* copyDeviceParas);
//double CalcuDeviceDist(Vector2 pos1, Vector2 pos2);
static __device__ int FindAxisIndex(double axis, const double* axisList, int axisCount);

//顺时针旋转后的坐标
static __device__ Vector2 Rotate(Vector2 pointPos, Vector2 centerPos, float rotateAngle);
//对一个数字*10000然后四舍五入到int
static __device__ int Multi10000ToInt(double num);

static __device__ Vector2Int Multi10000ToInt(Vector2 v);
//适应度计算函数 GPU
static __global__ void FitnessFunction(int curIterNum, int maxIterNum, int particleNum, int* bestParticleIndex,
	/*ProblemParas proParas, 固定参数的，不用管*/
	int DeviceSum, int fixedLinkPointSum, int fixedUniqueLinkPointSum, int vertPointCount, int horiPointCount, double workShopLength, double workShopWidth, double convey2DeviceDist, /*double conveyWidth, */
	double strConveyorUnitCost, double curveConveyorUnitCost, double conveyMinDist, /*double conveyMinLength, */double conveySpeed, Vector2 entrancePos, Vector2 exitPos,
	int CargoTypeNum, int totalLinkSum,

	/*CargoType* 固定参数*/
	/*int* deviceSum, */int* linkSum, int* accumLinkSum, DeviceLink* deviceLinkList, double* totalVolume,

	/*DevicePara**/
	Vector2* size, double* spaceLength, int* adjPInCount, int* adjPOutCount, int* accumAdjPInCount, int* accumAdjPOutCount,
	int totalInPoint, int totalOutPoint, AdjPoint* adjPointsIn, AdjPoint* adjPointsOut,
	/*Particle*/
	int dim, int fitnessCount, double* fitness_GPU, double* position_GPU, /*double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU*/
	/*存储所有粒子输送线路信息*/
	double* curBestFitnessVal, int inoutPSize, InoutPoint* inoutPoints, StraightConveyorInfo* strConveyorList,
	int* strConveyorListSum, Vector2Int* curveConveyorList, int* curveConveyorListSum,
	int* pointDirectArray, curandState* globalState)
{
	//粒子的下标i需要自己计算
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	double punishValue1 = 0;
	double punishValue2 = 0;
	bool IsDeviceOverlap = false;//是否重叠
	double deviceDist = 0;
	fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = 0;

#pragma region 深拷贝一份设备参数
	//数目：DeviceSum
	//DevicePara* copyDeviceParas_copy = new DevicePara[DeviceSum];
	//只有size需要copy
	Vector2* size_copy = new Vector2[DeviceSum];					//设备尺寸（分别是x轴和y轴的长度）

	//只有size需要复制一份新的
	for (int i = 0; i < DeviceSum; i++)
	{
		size_copy[i] = size[i];
	}
#pragma endregion

#pragma region 根据设备朝向，调整设备尺寸xy和出入口坐标
	int curAdjPIn_Index = 0;
	int curAdjPOut_Index = 0;
	for (int i = 2; i < dim; i += 3)
	{
		//double转int，转换为Direction，然后根据朝向重新计算设备尺寸和出入口
		//Rotate90或者Rotate270，尺寸的x和y互换
		//出入口按照顺时针算,旋转角=Direction*90(正好对应0,90,180,270）
		DeviceDirect curDirect = (DeviceDirect)(int)position_GPU[index * dim + i];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			double tempX = size_copy[i / 3].x;
			size_copy[i / 3].x = size_copy[i / 3].y;
			size_copy[i / 3].y = tempX;
		}
		//重新计算旋转后的出入口坐标
		Vector2 deviceCenterPos(0, 0, -1);//保证调用的是device版本的
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

#pragma region 检查设备是否重叠
	//如果重叠，进行调整
	//降低标准会发生什么？
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
				//cout << curIterNum << ":" << "重叠" << endl;
				//return;
			}
		}
	}
#pragma endregion

#pragma region 如果设备重叠，调整设备位置，尽力减少重叠，否则直接布线

#pragma region 调整设备位置
	if (IsDeviceOverlap == true) {
		//1.用一个设备的尺寸数组存储所有设备
		//2.每次从第一个设备开始，并且从编号为1的设备开始检测是否和该设备重叠
		//3.如果重叠，执行以下操作：
		//	a.规定设备只能向左下移动
		//  b.分别计算该设备向左/下移动的距离
		//	c.原则：选择移动的两个设备的尺寸x之和/尺寸y之和中相对于车间xy尺寸比例最小的方向移动
		//目前是移动到紧贴的位置，先不考虑溢出
		//4.移动之后，下次还是从第一个设备开始检测，直到所有的都不重叠
		//5.然后计算出所有设备位置的包络矩形，如果小于车间尺寸，移动到车间内的一个随机位置
		//然后根据这个移动前后的位置差修改所有设备的位置
		//用比例计算这个并不是特别靠谱，可以试试维护一个大的包络矩形
		int deviceIDSizeCount = DeviceSum;
		DeviceIDSize* deviceIDSizeList = new DeviceIDSize[DeviceSum];//按照设备大小排序的ID数组
		//先用其他结构存设备坐标，因为可能会修改失败
		double* particlePosList = new double[dim];
		for (int i = 0; i < dim; ++i) {
			particlePosList[i] = position_GPU[index * dim + i];
		}
		for (int i = 0; i < deviceIDSizeCount; ++i) {
			deviceIDSizeList[i] = DeviceIDSize(i, size_copy[i], -1);
		}
		DeviceIDSize_Sort(deviceIDSizeList, 0, deviceIDSizeCount - 1);//按照设备的尺寸排序

		double outSizeLength1, outSizeWidth1;
		double outSizeLength2, outSizeWidth2;
		int firstID, secondID;
		int maxIter = 1000;//防止死循环
		int curIter = 0;
		bool tooMuch = false;
		for (int i = 0; i < deviceIDSizeCount; ++i) {
			//检测其他的设备是否和它重叠
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
						j = 0;//只要遇到重叠的就要重新开始判断
						//重叠了，移动设备，根据比例确定往哪移动
						//这个策略不一定很好，需要验证
						double rateLeft = (outSizeLength1 + outSizeLength2) * 2 / workShopLength;
						double rateDown = (outSizeWidth1 + outSizeWidth2) * 2 / workShopWidth;
						if (rateLeft < rateDown) {//说明应该往左移动
							particlePosList[firstID * 3] = secondLowX - outSizeLength1 - 0.1;
						}
						else {//往下移动
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
		//5.然后计算出所有设备位置的包络矩形，如果小于车间尺寸，移动到车间内的一个随机位置
		//实现：计算出所有设备的四个方向的边界
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
		Vector2 oriRectAxis((max_X + min_X) / 2.0, (max_Y + min_Y) / 2.0, -1);//总包络矩形的中心坐标
		Vector2 newRectAxis;
		if ((max_X - min_X) <= workShopLength
			&& (max_Y - min_Y) <= workShopWidth) {//包络矩形小于车间大小
			//包络矩形随机产生一个坐标
			double rectLowX = 0 + (max_X - min_X) * 0.5;
			double rectHighX = workShopLength - (max_X - min_X) * 0.5;
			double rectLowY = 0 + (max_Y - min_Y) * 0.5;
			double rectHighY = workShopWidth - (max_Y - min_Y) * 0.5;
			newRectAxis.x = createARandomNum(globalState, index) * (rectHighX - rectLowX) + rectLowX;
			newRectAxis.y = createARandomNum(globalState, index) * (rectHighY - rectLowY) + rectLowY;
			//根据包络矩形的坐标变化，修改所有设备的坐标
			double deviceOffsetX = newRectAxis.x - oriRectAxis.x;
			double deviceOffsetY = newRectAxis.y - oriRectAxis.y;
			for (int i = 0; i < DeviceSum; ++i) {
				particlePosList[3 * i] += deviceOffsetX;
				particlePosList[3 * i + 1] += deviceOffsetY;
			}
			IsDeviceOverlap = false;
			//修改postion，需要修改其他的吗？
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
#pragma region 对齐方案1：对齐设备中心点x和y
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
		////对齐入口
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

#pragma region 对齐方案2：检测设备出入口点坐标并进行对齐操作
		//遍历所有cargoTypeList
		int curDeviceLink_Index = 0;//总的下标
		for (int j = 0; j < CargoTypeNum; j++)
		{
			for (int k = 0; k < linkSum[j]; k++)
			{
				int outDeviceIndex = deviceLinkList[curDeviceLink_Index].outDeviceIndex;
				int inDeviceIndex = deviceLinkList[curDeviceLink_Index].inDeviceIndex;
				int outPointIndex = deviceLinkList[curDeviceLink_Index].outPointIndex;
				int inPointIndex = deviceLinkList[curDeviceLink_Index].inPointIndex;
				AdjPoint outPoint, inPoint;
				//特殊情况1：仓库入口
				if (outDeviceIndex == -1)
				{
					//这里有个问题：下标的定位
					//每个adjPoint数组的大小不一样，这里可以设置一个记录前面x个adj数组点数目的数组
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (inPointTPos.x != entrancePos.x && abs(inPointTPos.x - entrancePos.x) < conveyMinDist)
					{
						//只能修改in，不能修改入口
						double moveLength = inPointTPos.x - entrancePos.x;
						position_GPU[index * dim + inDeviceIndex * 3] -= moveLength;

					}
					else if (inPointTPos.y != entrancePos.y && abs(inPointTPos.y - entrancePos.y) < conveyMinDist)
					{
						double moveLength = inPointTPos.y - entrancePos.y;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else if (inDeviceIndex == -2)//特殊情况2：仓库出口
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != exitPos.x && abs(outPointTPos.x - exitPos.x) < conveyMinDist)
					{
						//只能修改out，不能修改出口
						double moveLength = outPointTPos.x - exitPos.x;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
					}
					else if (outPointTPos.y != exitPos.y && abs(outPointTPos.y - exitPos.y) < conveyMinDist)
					{
						double moveLength = outPointTPos.y - exitPos.y;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else//其他情况
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					//在考虑了设备坐标的情况下对比
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
					{
						//x坐标接近
						double moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
						position_GPU[index * dim + inDeviceIndex * 3] += moveLength;

					}
					else if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)
					{
						//y坐标接近
						double moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] += moveLength;
					}
				}

				curDeviceLink_Index++;
			}
		}
#pragma endregion

#pragma region 对齐方案3：根据配对的出入口点的朝向分情况优化
		//遍历所有cargoTypeList
		//if 有一个是出口为i设备，且不是最后一个，那么就可以拿出这一对出入口点
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
				//特殊情况1：仓库入口
				if (outDeviceIndex == -1)
				{
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					if (inPointTPos.x != entrancePos.x && abs(inPointTPos.x - entrancePos.x) < conveyMinDist)
					{
						//只能修改in，不能修改入口
						double moveLength = inPointTPos.x - entrancePos.x;
						position_GPU[index * dim + inDeviceIndex * 3] -= moveLength;

					}
					else if (inPointTPos.y != entrancePos.y && abs(inPointTPos.y - entrancePos.y) < conveyMinDist)
					{
						double moveLength = inPointTPos.y - entrancePos.y;
						position_GPU[index * dim + inDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else if (inDeviceIndex == -2)//特殊情况2：仓库出口
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					if (outPointTPos.x != exitPos.x && abs(outPointTPos.x - exitPos.x) < conveyMinDist)
					{
						//只能修改out，不能修改出口
						double moveLength = outPointTPos.x - exitPos.x;
						position_GPU[index * dim + outDeviceIndex * 3] -= moveLength;
					}
					else if (outPointTPos.y != exitPos.y && abs(outPointTPos.y - exitPos.y) < conveyMinDist)
					{
						double moveLength = outPointTPos.y - exitPos.y;
						position_GPU[index * dim + outDeviceIndex * 3 + 1] -= moveLength;
					}
				}
				else//其他情况
				{
					outPoint = adjPointsOut[accumAdjPOutCount[outDeviceIndex] + outPointIndex];
					inPoint = adjPointsIn[accumAdjPInCount[inDeviceIndex] + inPointIndex];
					//考虑点的朝向
					outPointDirect = outPoint.direct;
					inPointDirect = inPoint.direct;
					//出入口的真实坐标
					Vector2 outPointTPos(outPoint.pos.x + position_GPU[index * dim + outDeviceIndex * 3],
						outPoint.pos.y + position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 inPointTPos(inPoint.pos.x + position_GPU[index * dim + inDeviceIndex * 3],
						inPoint.pos.y + position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					//出入口设备四个方向不考虑外层的边界
					Vector2 outDeviceUpPos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1] + 0.5 * size_copy[outDeviceIndex].y, -1);
					Vector2 outDeviceDownPos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1] - 0.5 * size_copy[outDeviceIndex].y, -1);
					Vector2 outDeviceLeftPos(position_GPU[index * dim + outDeviceIndex * 3] - 0.5 * size_copy[outDeviceIndex].x, position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					Vector2 outDeviceRightPos(position_GPU[index * dim + outDeviceIndex * 3] + 0.5 * size_copy[outDeviceIndex].x, position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);

					Vector2 inDeviceUpPos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1] + 0.5 * size_copy[inDeviceIndex].y, -1);
					Vector2 inDeviceDownPos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1] - 0.5 * size_copy[inDeviceIndex].y, -1);
					Vector2 inDeviceLeftPos(position_GPU[index * dim + inDeviceIndex * 3] - 0.5 * size_copy[inDeviceIndex].x, position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					Vector2 inDeviceRightPos(position_GPU[index * dim + inDeviceIndex * 3] + 0.5 * size_copy[inDeviceIndex].x, position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);


					//出入口设备四个方向的边界
					Vector2 outUpPos(outDeviceUpPos.x, outDeviceUpPos.y + convey2DeviceDist, -1);
					Vector2 outDownPos(outDeviceDownPos.x, outDeviceDownPos.y - convey2DeviceDist, -1);
					Vector2 outLeftPos(outDeviceLeftPos.x - convey2DeviceDist, outDeviceLeftPos.y, -1);
					Vector2 outRightPos(outDeviceRightPos.x + convey2DeviceDist, outDeviceRightPos.y, -1);

					Vector2 inUpPos(inDeviceUpPos.x, inDeviceUpPos.y + convey2DeviceDist, -1);
					Vector2 inDownPos(inDeviceDownPos.x, inDeviceDownPos.y - convey2DeviceDist, -1);
					Vector2 inLeftPos(inDeviceLeftPos.x - convey2DeviceDist, inDeviceLeftPos.x, -1);
					Vector2 inRightPos(outDeviceRightPos.x + convey2DeviceDist, outDeviceRightPos.y, -1);


					//要用到的各种比较值
					double inR2outXDist, out2inLXDist, in2outLXDist, outR2inXDist;
					double outU2inDYDist, inR2outLXDist, outR2inLXDist;
					double outL2inLXDist, inU2outUYDist, inU2outDYDist, outU2inYDist;
					double out2inDYDist, inU2outYDist, in2outDYDist;
					//先判断两者的方位
					//加一个判断二者方位的操作
					Vector2 inDevicePos(position_GPU[index * dim + inDeviceIndex * 3], position_GPU[index * dim + inDeviceIndex * 3 + 1], -1);
					Vector2 outDevicePos(position_GPU[index * dim + outDeviceIndex * 3], position_GPU[index * dim + outDeviceIndex * 3 + 1], -1);
					bool inLeft2OutBool, inRight2OutBool, inUp2OutBool, inDown2OutBool;

					int OutPosIndex_X, OutPosIndex_Y, InPosIndex_X, InPosIndex_Y;
					OutPosIndex_X = outDeviceIndex * 3; OutPosIndex_Y = outDeviceIndex * 3 + 1;
					InPosIndex_X = inDeviceIndex * 3; InPosIndex_Y = inDeviceIndex * 3 + 1;
					double moveLength = 0.0;
					switch (pointDirectArray[outPointDirect * 5 + inPointDirect])
					{
					case 1://上上
					{
						//1.两者在y上过于接近,让两者在y上对齐
						if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)
						{
							moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
							position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
							position_GPU[index * dim + InPosIndex_Y] += moveLength;
							break;
						}
						//2.入口设备在出口设备左/右上角
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
						//3.入口设备在出口设备左/右下角
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
					case 11://下下
					{
						//1.两者在y上过于接近,让两者在y上对齐
						if (outPointTPos.y != inPointTPos.y && abs(outPointTPos.y - inPointTPos.y) < conveyMinDist)//这个待定
						{
							moveLength = (outPointTPos.y - inPointTPos.y) * 0.5;
							position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
							position_GPU[index * dim + InPosIndex_Y] += moveLength;
							break;
						}
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						//2.入口设备在出口设备左/右上角
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
						//3.入口设备在出口设备左/右下角
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
					case 16://左左
					{
						//1.两者在x上过于接近,让两者在x上对齐
						if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
						{
							moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
							position_GPU[index * dim + OutPosIndex_X] -= moveLength;
							position_GPU[index * dim + InPosIndex_X] += moveLength;
							break;
						}
						//2.in在out左边
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//in在out左上
							out2inDYDist = outPointTPos.y - inDownPos.y;
							if (out2inDYDist > 0 && out2inDYDist < conveyMinDist)
							{
								moveLength = out2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in在out左下
							inU2outYDist = inUpPos.y - outPointTPos.y;
							if (inU2outYDist > 0 && inU2outYDist < conveyMinDist)
							{
								moveLength = inU2outYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						//3.in在out右边
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in在out右上
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in在out右下
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
					case 6://右右
					{
						//1.两者在x上过于接近,让两者在x上对齐
						if (outPointTPos.x != inPointTPos.x && abs(outPointTPos.x - inPointTPos.x) < conveyMinDist)
						{
							moveLength = (outPointTPos.x - inPointTPos.x) * 0.5;
							position_GPU[index * dim + OutPosIndex_X] -= moveLength;
							position_GPU[index * dim + InPosIndex_X] += moveLength;
							break;
						}
						//2.in在out右边
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in在out右上
							out2inDYDist = outPointTPos.y - inDownPos.y;
							if (out2inDYDist > 0 && out2inDYDist < conveyMinDist)
							{
								moveLength = out2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in在out右下
							inU2outYDist = inUpPos.y - outPointTPos.y;
							if (inU2outYDist > 0 && inU2outYDist < conveyMinDist)
							{
								moveLength = inU2outYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] += moveLength;
								position_GPU[index * dim + InPosIndex_Y] -= moveLength;
								break;
							}
						}
						//3.in在out左边
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//in在out左上
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
							//in在out左下
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
					case 3://上下
					{
						//1.in在out上方，且在x上很接近
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
						//2.入口设备在出口设备左/右上角
						inUp2OutBool = IsInUp2Out(inDeviceDownPos, outDeviceUpPos);
						if (inUp2OutBool)//左右可以用一种方式计算
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
						//入口设备在出口设备下面
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							//3.入口设备在出口设备下面/分为左右
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
					case 14://左右
					{
						//1.in在out左边，且在y上很接近
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
						//2.in在out左边
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
						//3.in在out右边，分为右上和右下
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);
						if (inRight2OutBool)
						{
							//in在out右下
							inU2outDYDist = inUpPos.y - outDownPos.y;
							if (inU2outDYDist > 0 && inU2outDYDist < conveyMinDist)
							{
								moveLength = inU2outDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
							//in在out右上
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
					case 4://上左(先不考虑可能的情况)
					{
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//1.in在out左上
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
							//2.in在out右上
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
							//3.in在out右上一点
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
							//4.in在out右下角
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
					case 2://上右
					{
						//in在out上
						inUp2OutBool = IsInUp2Out(inDownPos, outUpPos);
						if (inUp2OutBool)
						{
							//1.in在out左上
							inR2outXDist = inRightPos.x - outPointTPos.x;
							if (inR2outXDist > 0 && inR2outXDist < conveyMinDist)
							{
								moveLength = inR2outXDist * 0.5;
								position_GPU[index * dim + OutPosIndex_X] += moveLength;
								position_GPU[index * dim + InPosIndex_X] -= moveLength;
								break;
							}
						}
						//in在out右
						inRight2OutBool = IsInRight2Out(inLeftPos, outRightPos);//这个计算也可以优化，提前算好
						if (inRight2OutBool)
						{
							//2.in在out右上
							outU2inDYDist = outUpPos.y - inDownPos.y;
							if (outU2inDYDist > 0 && outU2inDYDist < conveyMinDist)
							{
								moveLength = outU2inDYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//in在out左
						inLeft2OutBool = IsInLeft2Out(inRightPos, outLeftPos);
						if (inLeft2OutBool)
						{
							//3.in在out左上
							outU2inYDist = outUpPos.y - inPointTPos.y;
							if (outU2inYDist > 0 && outU2inYDist < conveyMinDist)
							{
								moveLength = outU2inYDist * 0.5;
								position_GPU[index * dim + OutPosIndex_Y] -= moveLength;
								position_GPU[index * dim + InPosIndex_Y] += moveLength;
								break;
							}
						}
						//in在out下
						inDown2OutBool = IsInDown2Out(inUpPos, outDownPos);
						if (inDown2OutBool)
						{
							//3.in在out左下
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

#pragma region 计算出入口点的集合坐标
		//现在in和out点被分开了
		//然后给所有的inoutPoint赋值
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
				inoutPoints[index * inoutPSize + ioPIndex] = ioPoint;//注意偏移值
				ioPIndex++;
			}
			for (int pointIndex = 0; pointIndex < adjPOutCount[i]; ++pointIndex)
			{
				AdjPoint& point = adjPointsOut[accumAdjPOutCount[i] + pointIndex];
				InoutPoint ioPoint;
				ioPoint.pointDirect = point.direct;
				Vector2 axis(point.pos.x + position_GPU[index * dim + 3 * i], point.pos.y + position_GPU[index * dim + 3 * i + 1], -1);
				ioPoint.pointAxis = axis;
				inoutPoints[index * inoutPSize + ioPIndex] = ioPoint;//注意偏移值
				ioPIndex++;
			}
		}

#pragma endregion

#pragma region 根据设备坐标和出入口坐标构造路径点图
		double* horizonAxisList = new double[horiPointCount];
		double* verticalAxisList = new double[vertPointCount];
		int curHoriIndex = 0;
		int curVertIndex = 0;
		//先对出入口点水平和垂直进行分类(注意加上偏移量)
		//先求出水平和垂直出入口点的数目（一部分horiCount和vertCount）
		for (int i = 0; i < DeviceSum; i++)
		{
			for (int pointIndex = 0; pointIndex < adjPInCount[i]; pointIndex++)
			{
				AdjPoint& p = adjPointsIn[accumAdjPInCount[i] + pointIndex];
				if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//上下
				{
					horizonAxisList[curHoriIndex++] = p.pos.x + position_GPU[index * dim + i * 3];
				}
				else {//左右
					verticalAxisList[curVertIndex++] = p.pos.y + position_GPU[index * dim + i * 3 + 1];
				}
			}
			for (int pointIndex = 0; pointIndex < adjPOutCount[i]; pointIndex++)
			{
				AdjPoint& p = adjPointsOut[accumAdjPOutCount[i] + pointIndex];
				if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//上下
				{
					horizonAxisList[curHoriIndex++] = p.pos.x + position_GPU[index * dim + i * 3];
				}
				else {//左右
					verticalAxisList[curVertIndex++] = p.pos.y + position_GPU[index * dim + i * 3 + 1];
				}
			}
		}
		//仓库入口2个点&出口（horiCount和vertCount+=2）
		horizonAxisList[curHoriIndex++] = entrancePos.x;
		verticalAxisList[curVertIndex++] = entrancePos.y;

		horizonAxisList[curHoriIndex++] = exitPos.x;
		verticalAxisList[curVertIndex++] = exitPos.y;
		//存下每个设备坐标的四个范围（作为后面障碍点的范围）
		double* DeviceLowXList = new double[DeviceSum];
		double* DeviceHighXList = new double[DeviceSum];
		double* DeviceLowYList = new double[DeviceSum];
		double* DeviceHighYList = new double[DeviceSum];
		//每个设备周围的4个点
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

			//每个设备的四个范围
			DeviceLowXList[i / 3] = LowX;
			DeviceHighXList[i / 3] = HighX;
			DeviceLowYList[i / 3] = LowY;
			DeviceHighYList[i / 3] = HighY;

			//防止路径进入设备内部
			horizonAxisList[curHoriIndex++] = position_GPU[index * dim + i];
			verticalAxisList[curVertIndex++] = position_GPU[index * dim + i + 1];

		}
		//进一步处理这些点
		//对这些点的坐标按照从小到大排序
		//sort函数需要自己实现，或者使用cuda库函数的
		Double_Sort(horizonAxisList, 0, horiPointCount - 1);
		Double_Sort(verticalAxisList, 0, vertPointCount - 1);
		//只保留不重复的点（自己实现）

		int uniqueHoriPCount = horiPointCount;
		int uniqueVertPCount = vertPointCount;
		int unique_end1 = Double_Unique(horizonAxisList, 0, uniqueHoriPCount - 1);
		uniqueHoriPCount = unique_end1;
		//verticalAxisList.erase(unique_end1, verticalAxisList.end());
		int unique_end2 = Double_Unique(verticalAxisList, 0, uniqueVertPCount - 1);
		uniqueVertPCount = unique_end2;
		//horizonAxisList.erase(unique_end2, horizonAxisList.end());

		//存所有的障碍点的下标
		int barrierRowNum = 200;//用一个固定大小分配内存
		int* barrierRowIndexList = new int[barrierRowNum];
		int barrierColNum = 200;
		int* barrierColIndexList = new int[barrierColNum];
		int totalBarRowNum = 0;//记录barrier的实际行数目
		int totalBarColNum = 0;//记录barrier的实际列数目

		//用这些坐标去组成路径点map，map是二维的，相当于二维点矩阵
		//用一维代替二维
		//horiNum对应colNum，vertNum对应rowNum
		int pathColNum = uniqueHoriPCount;
		int pathRowNum = uniqueVertPCount;
		APoint** pathPointMap = new APoint*[pathColNum * pathRowNum];

		for (int i = 0; i < pathColNum * pathRowNum; i++)
		{
			pathPointMap[i] = new APoint();//二级分配内存
			int rowIndex = i / pathColNum;//对应下面的i
			int colIndex = i % pathColNum;//对应下面的j
			pathPointMap[i]->x = horizonAxisList[colIndex];
			pathPointMap[i]->y = verticalAxisList[rowIndex];
			//遍历所有设备，看是否有和这个点重叠的（实现标记障碍点）
			for (int k = 0; k < DeviceSum; k++)
			{
				if (pathPointMap[i]->x - DeviceLowXList[k] >= 0.01 && DeviceHighXList[k] - pathPointMap[i]->x >= 0.01
					&& pathPointMap[i]->y - DeviceLowYList[k] >= 0.01 && DeviceHighYList[k] - pathPointMap[i]->y >= 0.01)
				{
					pathPointMap[i]->type = AType::ATYPE_BARRIER;
					//障碍点在图中的下标
					barrierRowIndexList[totalBarRowNum++] = rowIndex;
					barrierColIndexList[totalBarColNum++] = colIndex;
				}
			}
			pathPointMap[i]->colIndex = colIndex;//记录下点在路径点map中的下标
			pathPointMap[i]->rowIndex = rowIndex;
		}
#pragma endregion

#pragma region 寻路
		CAstar* star = new CAstar();
		//只在kernal函数内部执行，用new和delete就行
		star->_allPoints = pathPointMap;
		star->pointColNum = pathColNum;
		star->pointRowNum = pathRowNum;

		int beginRowIndex, beginColIndex, endRowIndex, endColIndex;

		double totalTime = 0.0;
		//存所有的路径信息
		PointLink* copyPLinks = new PointLink[totalLinkSum];//此时每条路径有50个点
		int curLinkIndex = 0;//当前的link下标
		int totalLinkPointSum = 0;//所有link的所有点的数目（每个link有很多点）
		for (int i = 0; i < CargoTypeNum; i++)//货物类型
		{
			//CargoType curCargoType = proParas.cargoTypeList[i];
			//每种货物可能经过多个设备，也就是一个Type对应多个link
			for (int j = 0; j < linkSum[i]; j++)
			{
				PathDirection pathBeginDirect;
				PathDirection pathEndDirect;
				int forwardDeviceIndex, curDeviceIndex;//设备1和设备2ID
				int forwardOutIndex, curInIndex;//出入口的下标
				double device1PosX, device1PosY, device2PosX, device2PosY;//设备周围的四个点
				double initDevice1PosX, initDevice1PosY, initDevice2PosX, initDevice2PosY;//保存未增加包围边的坐标

				double deviceDistance = 0.0;//距离
				double outDSizeL, inDSizeL;

				forwardDeviceIndex = deviceLinkList[accumLinkSum[i] + j].outDeviceIndex;
				curDeviceIndex = deviceLinkList[accumLinkSum[i] + j].inDeviceIndex;

				forwardOutIndex = deviceLinkList[accumLinkSum[i] + j].outPointIndex;
				curInIndex = deviceLinkList[accumLinkSum[i] + j].inPointIndex;
				if (forwardDeviceIndex == -1)//说明是仓库入口
				{
					/*forwardOutIndex = 0;
					curInIndex = proParas.cargoTypeList[i].deviceLinkList[j].inPointIndex;*/

					//开头和结尾点的朝向方向
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
					//得到设备周围的点
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
				else if (curDeviceIndex == -2)//说明是出口
				{
					//开头和结尾点的朝向方向 
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
					//得到设备周围的点 
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
				else//普通
				{
					//forwardOutIndex = proParas.cargoTypeList[i].deviceLinkList[j].outPointIndex - 1;
					//curInIndex = proParas.cargoTypeList[i].deviceLinkList[j].inPointIndex - 1;

					//开头和结尾点的朝向方向
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
					//得到设备周围的点 
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
				//计算最短路径
				beginRowIndex = FindAxisIndex(device1PosY, verticalAxisList, uniqueVertPCount);
				beginColIndex = FindAxisIndex(device1PosX, horizonAxisList, uniqueHoriPCount);
				endRowIndex = FindAxisIndex(device2PosY, verticalAxisList, uniqueVertPCount);
				endColIndex = FindAxisIndex(device2PosX, horizonAxisList, uniqueHoriPCount);

				//得到路径，path是第一个节点
				APoint* path = findWay(star->curPathDirect, pathBeginDirect, star->_allPoints, star->_endPoint, star->_neighbourList,
					star->_curPoint, star->_openList, star->_closeList, star->openList_CurSize, star->closeList_CurSize, star->neighbourList_CurSize,
					star->pointColNum, star->pointRowNum, beginRowIndex, beginColIndex, endRowIndex, endColIndex);
				//不可行的解，直接退出
				if (path == nullptr)
				{
					fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
					printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
					return;
				}
				//根据路径计算长度
				deviceDistance = CalcuPathLength(path) + outDSizeL + inDSizeL;


#pragma region 计算路径（只带上起始点 + 路径中的转弯点）
				//路径保存下来
				//这个路径中的点还没有简化，先给每个路径设置大小为50
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
				points1[points1Index++] = startP1;//points1Index是某一个link的实际点的数目
				totalLinkPointSum += points1Index;//统计link中所有的实际点的数目（方便之后分配空间）

				PointLink pointLink1(forwardDeviceIndex, forwardOutIndex, curDeviceIndex, curInIndex, points1, points1Index, -1);
				copyPLinks[curLinkIndex++] = pointLink1;

				//Vector2 lastP(path->x, path->y);//从最后一个点开始
				//points.push_back(lastP);
				//path = path->parent;
				//while (path)
				//{
				//	Vector2 curP(path->x, path->y);
				//	//只有转弯的点才会被加入路径
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

				//计算输送时间(物料总量 * 路线长度 * 输送效率) 
				totalTime += totalVolume[i] * deviceDistance * conveySpeed;
				//计算设备处理时间(物料总量 * 处理效率)
				//totalTime += curCargoType.totalVolume * curDevice.workSpeed;

				resetAStar(star->pointRowNum, star->pointColNum, star->_allPoints, star->_openList, star->_closeList, star->_neighbourList,
					star->openList_CurSize, star->closeList_CurSize, star->neighbourList_CurSize);
				//给障碍点重新标记
				for (int i = 0; i < totalBarRowNum; i++)
				{
					star->_allPoints[barrierRowIndexList[i] * pathColNum + barrierColIndexList[i]]->type = AType::ATYPE_BARRIER;
				}
			}
		}
#pragma endregion

#pragma region 将布局结果转化为输送机参数(同时加强一下输送线最短距离约束)

		//只保留转弯点的数组
		int segPathSet_PointSum = fixedLinkPointSum * totalLinkSum;//50个点
		int segPathSet_CurIndex = 0;
		SegPath* segPathSet = new SegPath[segPathSet_PointSum];

		//////过滤所有重复的路线段，获取segPathSet//////
		//1.将copyPLinks里的所有点放到seg数组里
		for (int i = 0; i < totalLinkSum; i++)
		{
			for (int j = copyPLinks[i].pointNum - 1; j > 0; j--)
			{
				Vector2Int p1 = Multi10000ToInt(copyPLinks[i].points[j]);
				Vector2Int p2 = Multi10000ToInt(copyPLinks[i].points[j - 1]);
				if (p1.ANEqualB(p2, -1))//坐标不能一样
				{
					SegPath temp(p1, p2, -1);
					segPathSet[segPathSet_CurIndex++] = temp;
				}
			}
		}
		//2.对seg数组进行排序+去重
		SegPath_Sort(segPathSet, 0, segPathSet_CurIndex - 1);
		//去重后的点的数目
		int segPathSet_UniqueSum = SegPath_Unique(segPathSet, 0, segPathSet_CurIndex - 1);

		//////遍历set，计算出set中所有的 点所在路线的垂直水平数目（类似出入度）&是否保留 信息//////
		//这些信息存到pathPointInfoMap中
		//map<Vector2Int, PointInfo> pathPointInfoMap;


		//路径点信息map(每个点的坐标:每个点的信息)
		int pathPointInfoMap_PointSum = segPathSet_UniqueSum * 2;//每个seg对应两个点
		PointInfo* pathPointInfoMap = new PointInfo[pathPointInfoMap_PointSum];
		int pathPointInfoMap_CurIndex = 0;//当前下标

		//同样的，只能用数组代替map
		//如何实现？
		//1.现将segPath里所有的点放到PointInfo里
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
		//2.对所有点进行排序，目的是将所有相同的点聚集到一起
		PointInfo_Sort(pathPointInfoMap, 0, pathPointInfoMap_CurIndex - 1);
		//3.根据排序的数组，更新每个点的vertNum和horiNum数目,顺便去重
		int pathPointInfoMap_UniqueSum = PointInfo_CalcuAndUnique(pathPointInfoMap, 0, pathPointInfoMap_CurIndex - 1);


		//4.根据点的vertNum和horiNum，判断每个点是否被保留
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

		//从pathPointInfoMap中得到tempStrConveyorList和tempCurveConveyorList
		//原来由于有set和map，可以实现快读查找
		//现在没有了，只能用二分查找了

		//直线输送机信息列表
		int tempStrConveyorList_PointSum = fixedUniqueLinkPointSum * totalLinkSum;//20个点
		//StraightConveyorInfo* tempStrConveyorList = new StraightConveyorInfo[tempStrConveyorList_PointSum];
		int tempStrConveyorList_CurIndex = 0;

		//转弯输送机信息列表
		int tempCurveConveyorList_PointSum = fixedUniqueLinkPointSum * totalLinkSum;//20个点
		//Vector2Int* tempCurveConveyorList = new Vector2Int[tempCurveConveyorList_PointSum];
		int tempCurveConveyorList_CurIndex = 0;

		for (int i = 0; i < totalLinkSum; i++) {
			StraightConveyorInfo tempStrInfo;
			tempStrInfo.startPos = Multi10000ToInt(copyPLinks[i].points[copyPLinks[i].pointNum - 1]);//开头
			PointInfo curPointInfo = FindPointInfo(pathPointInfoMap, 0, pathPointInfoMap_UniqueSum - 1, tempStrInfo.startPos);
			tempStrInfo.startVnum = curPointInfo.vertDirNum;
			tempStrInfo.startHnum = curPointInfo.horiDirNum;
			for (int j = copyPLinks[i].pointNum - 2; j >= 0; j--) {
				Vector2Int p = Multi10000ToInt(copyPLinks[i].points[j]);
				PointInfo pPathPoint = FindPointInfo(pathPointInfoMap, 0, pathPointInfoMap_UniqueSum - 1, p);
				if (pPathPoint.isKeep == true) {//这里会出现重复的
					//先更新直线输送机
					if (tempStrInfo.startPos.ANEqualB(p, -1)) {
						tempStrInfo.endPos = p;
						//if (curIterNum != 0)
						if (tempStrInfo.startPos.Distance(tempStrInfo.endPos) < conveyMinDist)
						{
							//cout << "输送线太短" << endl;
							//punishValue1 = 150 * (curIterNum + 1);
							//punishValue2 = 30 * (curIterNum + 1);
							fitness_GPU[index * fitnessCount + 0] = fitness_GPU[index * fitnessCount + 1] = MAX_FITNESS;
							printf("%d\n", fitness_GPU[index * fitnessCount + 0]);
							return;
						}
						tempStrInfo.endVnum = pPathPoint.vertDirNum;
						tempStrInfo.endHnum = pPathPoint.horiDirNum;
						//这个地方进行去重的
						//如果没有set，只能先排序+去重复了
						//注意每个p既是终点也是下一个的起点
						strConveyorList[index * tempStrConveyorList_PointSum + tempStrConveyorList_CurIndex++] = tempStrInfo;//注意偏移值
						tempStrInfo.startPos = p;
						tempStrInfo.startVnum = pPathPoint.vertDirNum;
						tempStrInfo.startHnum = pPathPoint.horiDirNum;
					}
					//只要不是始终点，都需要更新转弯输送机
					if (!(pPathPoint.horiDirNum == 1 && pPathPoint.vertDirNum == 0)
						&& !(pPathPoint.horiDirNum == 0 && pPathPoint.vertDirNum == 1)) {
						//tempCurveConveyorList.insert(p);
						curveConveyorList[index * tempCurveConveyorList_PointSum + tempCurveConveyorList_CurIndex++] = p;
					}
				}
			}
		}
#pragma endregion

#pragma region 计算目标函数值
		//遍历直线和转弯输送机的set，得到输送机的总成本
		double conveyorTotalCost = 0.0;
		for (int i = 0; i < tempStrConveyorList_CurIndex; i++) {
			conveyorTotalCost += strConveyorUnitCost * strConveyorList[index * tempStrConveyorList_PointSum + i]
				.startPos.Distance(strConveyorList[index * tempStrConveyorList_PointSum + i].endPos);
		}
		conveyorTotalCost += curveConveyorUnitCost * tempStrConveyorList_CurIndex;
		//设置适应度值
		fitness_GPU[index * fitnessCount + 0] = totalTime;
		fitness_GPU[index * fitnessCount + 1] = conveyorTotalCost;

		//上面相当于计算出了可能的最佳输送线，现在需要更新bestPathInfoList
		//根据适应度是否升级选择更新BestPathInfoList
		//一个可以并行计算最大值的归约算法
		const int ParticleNum = 100;
		__shared__ int particleIndexList[ParticleNum];
		//为其赋值
		particleIndexList[index] = index;
		//等待每个线程赋值完成
		__syncthreads();
		//实现归约，查找所有粒子中最佳粒子的index(最佳index存到bestParticleIndex[0])
		int leng = particleNum;
		for (int i = particleNum / 2.0 + 0.5; i > 1; i = i / 2.0 + 0.5)
		{
			if (index < i)
			{
				if (index + i < leng)
				{
					//默认处理适应度1
					if (fitness_GPU[particleIndexList[index] * fitnessCount + 0]
						< fitness_GPU[(particleIndexList[index] + i) * fitnessCount + 0]) {
						//交换
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
		//bestParticleIndex[0]存的是最佳的fitness对应的粒子下标
		//for (int i = 0; i < fitnessCount; ++i) {
		//	if (fitness_GPU[index * fitnessCount + i] < bestPathInfoList[i].curBestFitnessVal) {//需要更新
		//		bestPathInfoList[i].curBestFitnessVal = fitness_GPU[index * fitnessCount + i];
		//		bestPathInfoList[i].inoutPoints = inoutPoints + index * inoutPSize;//注意偏移值
		//		bestPathInfoList[i].inoutPSize = totalInPoint + totalOutPoint;//可能有问题
		//		bestPathInfoList[i].strConveyorList = strConveyorList + index * tempStrConveyorList_PointSum;
		//		bestPathInfoList[i].strConveyorListSum = tempStrConveyorList_CurIndex;//大小要记录，初始分配大小：fixedUniqueLinkPointSum * totalLinkSum * fitnessCount
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
//顺时针旋转后的坐标
static __device__ Vector2 Rotate(Vector2 pointPos, Vector2 centerPos, float rotateAngle)
{
	float xx = (pointPos.x - centerPos.x) * cos(rotateAngle * (PI / 180)) + (pointPos.y - centerPos.y) * sin(rotateAngle * (PI / 180)) + centerPos.x;
	float yy = -(pointPos.x - centerPos.x) * sin(rotateAngle * (PI / 180)) + (pointPos.y - centerPos.y) * cos(rotateAngle * (PI / 180)) + centerPos.y;
	Vector2 result(xx, yy, -1);
	return result;
}
static __device__ int FindAxisIndex(double axis, const double* axisList, int axisCount)
{
	//用二分法更快
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
////计算占地面积
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
//	//计算总面积
//	area = (max_X - min_X) * (max_Y - min_Y);
//	return area;
//}
//先乘以10000，然后四舍五入到Int
static __device__ int Multi10000ToInt(double num)
{
	//使用自定义的Round函数
	return MyRound(num * 10000);
}
static __device__ Vector2Int Multi10000ToInt(Vector2 v)
{
	return Vector2Int(Multi10000ToInt(v.x), Multi10000ToInt(v.y), -1);
}