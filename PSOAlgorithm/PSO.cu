#pragma once
#include "Pareto.h"
#include "Archive.h"
#include "FitnessFunction.h"
#include "Tools.h"
#include <ctime>


//��дstatic

static __global__ void InitPointDirectArray(int* pointDirectArray)
{
	pointDirectArray[0] = -1;
	pointDirectArray[1] = -1;
	pointDirectArray[2] = -1;
	pointDirectArray[3] = -1;
	pointDirectArray[4] = -1;

	pointDirectArray[5] = -1;
	pointDirectArray[6] = 1;
	pointDirectArray[7] = 2;
	pointDirectArray[8] = 3;
	pointDirectArray[9] = 4;

	pointDirectArray[10] = -1;
	pointDirectArray[11] = 5;
	pointDirectArray[12] = 6;
	pointDirectArray[13] = 7;
	pointDirectArray[14] = 8;

	pointDirectArray[15] = -1;
	pointDirectArray[16] = 9;
	pointDirectArray[17] = 10;
	pointDirectArray[18] = 11;
	pointDirectArray[19] = 12;

	pointDirectArray[20] = -1;
	pointDirectArray[21] = 13;
	pointDirectArray[22] = 14;
	pointDirectArray[23] = 15;
	pointDirectArray[24] = 16;
}
static __global__ void InitCurBestPathFit(double* curBestPath_FitnessVal)
{
	curBestPath_FitnessVal[0] = 1000000;
}

//��������&������Ӧ��
static __global__ void UpdateParticle_Kernal(int curIterNum, int maxIterNum, int dim, int fitnessCount, double w_, double C1_, double C2_, double dt_,
	/*�������Particle Particle* particles_,*/double* fitness_GPU, double* position_GPU, double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU,
	curandState* globalState, double* randomNumList, double* range_interval_, double* upper_bound_, double* lower_bound_, double* all_best_position_,
	/*ProblemParas problemParas*/Vector2* size, double* spaceLength, double workShopLength, double workShopWidth)
{
	//���ӵ��±�i��Ҫ�Լ�����
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	//�ȸ��³���Ȼ����ݳ���������ӵķ�Χ
	for (int j = 2; j < dim; j += 3)
	{
		double last_position = position_GPU[i * dim + j];

		velocity_GPU[i * dim + j] = w_ * velocity_GPU[i * dim + j] +
			C1_ * createARandomNum(globalState, i) * (best_position_GPU[i * dim + j] - position_GPU[i * dim + j]) +
			C2_ * createARandomNum(globalState, i) * (all_best_position_[i * dim + j] - position_GPU[i * dim + j]);
		position_GPU[i * dim + j] += dt_ * velocity_GPU[i * dim + j];

		// �����������������������
		if (upper_bound_ && lower_bound_)
		{
			if (position_GPU[i * dim + j] >= upper_bound_[j])//ע������豸����=Ҳ����
			{
				double thre = createARandomNum(globalState, i);//ֱ������һ�������
				if (last_position >= upper_bound_[j] - 1)//ע��upper_bound_[j]-1=3
				{
					position_GPU[i * dim + j] = createARandomNum(globalState, i) * range_interval_[j] + lower_bound_[j];
				}
				else if (thre < 0.5)
				{
					position_GPU[i * dim + j] = upper_bound_[j] - (upper_bound_[j] - last_position) * createARandomNum(globalState, i);
				}
				else
				{
					position_GPU[i * dim + j] = upper_bound_[j] - 0.5;
				}
			}
			if (position_GPU[i * dim + j] < lower_bound_[j])
			{
				double thre = createARandomNum(globalState, i);
				if (last_position == lower_bound_[j])
				{
					position_GPU[i * dim + j] = createARandomNum(globalState, i) * range_interval_[j] + lower_bound_[j];
				}
				else if (thre < 0.5)
				{
					position_GPU[i * dim + j] = lower_bound_[j] + (last_position - lower_bound_[j]) * createARandomNum(globalState, i);
				}
				else
				{
					position_GPU[i * dim + j] = lower_bound_[j];
				}
			}
		}
	}
	//���ݳ����޸��豸���½緶Χ
	for (int j = 2; j < dim; j += 3)
	{
		//doubleתint��ת��ΪDirection��Ȼ����ݳ������¼����豸�ߴ�ͳ����
		//Rotate90����Rotate270,�޸�������
		DeviceDirect curDirect = (DeviceDirect)(int)position_GPU[i * dim + j];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)//��һ���ֿ���ҲҪ�ģ�enum��C++�﷨��
		{
			//x��y
			lower_bound_[j - 2] = 0 + size[j / 3].y * 0.5 + spaceLength[j / 3];
			lower_bound_[j - 1] = 0 + size[j / 3].x * 0.5 + spaceLength[j / 3];

			upper_bound_[j - 2] = workShopLength - size[j / 3].y * 0.5 - spaceLength[j / 3];
			upper_bound_[j - 1] = workShopWidth - size[j / 3].x * 0.5 - spaceLength[j / 3];

		}
		else
		{
			//x��y
			lower_bound_[j - 2] = 0 + size[j / 3].x * 0.5 + spaceLength[j / 3];
			lower_bound_[j - 1] = 0 + size[j / 3].y * 0.5 + spaceLength[j / 3];

			upper_bound_[j - 2] = workShopLength - size[j / 3].x * 0.5 - spaceLength[j / 3];
			upper_bound_[j - 1] = workShopWidth - size[j / 3].y * 0.5 - spaceLength[j / 3];

		}
		range_interval_[j - 2] = upper_bound_[j - 2] - lower_bound_[j - 2];
		range_interval_[j - 1] = upper_bound_[j - 1] - lower_bound_[j - 1];
	}
	//cout << endl;
	for (int j = 0; j < dim; j++)
	{
		if (j % 3 != 2)
		{
			//������һ�ε��������position��velocity
			double last_position = position_GPU[i * dim + j];

			velocity_GPU[i * dim + j] = w_ * velocity_GPU[i * dim + j] +
				C1_ * createARandomNum(globalState, i) * (best_position_GPU[i * dim + j] - position_GPU[i * dim + j]) +
				C2_ * createARandomNum(globalState, i) * (all_best_position_[i * dim + j] - position_GPU[i * dim + j]);
			position_GPU[i * dim + j] += dt_ * velocity_GPU[i * dim + j];

			// �����������������������
			if (upper_bound_ && lower_bound_)
			{
				if (position_GPU[i * dim + j] > upper_bound_[j])
				{
					double thre = createARandomNum(globalState, i);
					if (last_position >= upper_bound_[j])
					{
						position_GPU[i * dim + j] = createARandomNum(globalState, i) * range_interval_[j] + lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						position_GPU[i * dim + j] = upper_bound_[j] - abs(upper_bound_[j] - last_position) * createARandomNum(globalState, i);
					}
					else
					{
						position_GPU[i * dim + j] = upper_bound_[j];
					}
				}
				if (position_GPU[i * dim + j] < lower_bound_[j])
				{
					double thre = createARandomNum(globalState, i);
					if (last_position <= lower_bound_[j])
					{
						position_GPU[i * dim + j] = createARandomNum(globalState, i) * range_interval_[j] + lower_bound_[j];
					}
					else if (thre < 0.5)
					{
						position_GPU[i * dim + j] = lower_bound_[j] + abs(last_position - lower_bound_[j]) * createARandomNum(globalState, i);
					}
					else
					{
						position_GPU[i * dim + j] = lower_bound_[j];
					}
				}
			}
		}

	}
}
// �Ƚ��������ӵ���Ӧ�ȣ��ж��Ƿ���ȫ֧�䣬�Ӷ������pbest
static __device__ bool ComparePbest(int index, int fitness_count, double* fitness, double* pbestFitness, curandState* globalState)
{
	int numGreater = 0;
	int numLess = 0;
	for (int i = 0; i < fitness_count; i++)
	{
		if (fitness[i] < pbestFitness[i])
		{
			numGreater++;
		}
		if (fitness[i] > pbestFitness[i])
		{
			numLess++;
		}
	}
	//�����ǰ֧����ʷ������
	if (numGreater > 0 && numLess == 0)
	{
		return true;
	}
	//�����ʷ֧�䵱ǰ���ӣ�������
	else if (numGreater == 0 && numLess > 0)
	{
		return false;
	}
	//�������֧�䣬���ѡ����Ӧ��1�ĸ��ʸߵ㣩
	else
	{
		double randomProb = createARandomNum(globalState, index);//�������С��
		if (fitness[0] < pbestFitness[0])
		{
			if (randomProb < 0.75)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			if (randomProb > 0.5)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}

//����Pbest��GPU����
static __global__ void UpdatePbest_Kernal(int dim_, int fitness_count, double* fitness_GPU, double* position_GPU,
	double* velocity_GPU, double* best_position_GPU, double* best_fitness_GPU, curandState* globalState)
{
	//i��Ҫ����
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	//�Ƚ���ʷpbest�͵�ǰ��Ӧ�ȣ������Ƿ�Ҫ����
	if (ComparePbest(i, fitness_count, fitness_GPU + i * fitness_count, best_fitness_GPU + i * fitness_count, globalState));
	{
		for (int j = 0; j < fitness_count; j++)
		{
			best_fitness_GPU[i * fitness_count + j] = fitness_GPU[i * fitness_count + j];
			//particles_[i].best_fitness_[j] = particles_[i].fitness_[j];
		}
		for (int j = 0; j < dim_; j++)
		{
			best_position_GPU[i * dim_ + j] = position_GPU[i * dim_ + j];
			//particles_[i].best_position_[j] = particles_[i].position_[j];
		}
	}
}

//����Gbest��GPU����
static __global__ void UpdateGbest_Kernal(int fitness_count, int dim_, double* all_best_fitness_, double* all_best_position_, Particle* gbestList)
{
	//�±��Լ���
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	//����
	for (int j = 0; j < fitness_count; j++)
	{
		all_best_fitness_[i * fitness_count + j] = gbestList[i].best_fitness_[j];
	}
	for (int k = 0; k < dim_; k++)
	{
		all_best_position_[i * dim_ + k] = gbestList[i].best_position_[k];
	}
}





// ���캯��(��ʼ�������㷨�Ĳ��������������ռ�)
PSOOptimizer::PSOOptimizer(PSOPara* pso_para, ProblemParas& problemParas)
{
	blockNum = pso_para->blockSum;
	threadsPerBlock = pso_para->threadsPerBlock;

	particle_num_ = pso_para->particle_num_;
	max_iter_num_ = pso_para->max_iter_num_;
	dim_ = pso_para->dim_;
	fitness_count = pso_para->fitness_count_;
	curr_iter_ = 0;
	meshDivCount = pso_para->mesh_div_count;
	archiveMaxCount = pso_para->archive_max_count;

	//���Ӳ���CPU 
	position_CPU = (double*)malloc(sizeof(double) * particle_num_ * dim_);
	velocity_CPU = (double*)malloc(sizeof(double) * particle_num_ * dim_);
	best_position_CPU = (double*)malloc(sizeof(double) * particle_num_ * dim_);
	fitness_CPU = (double*)malloc(sizeof(double) * particle_num_ * fitness_count);
	best_fitness_CPU = (double*)malloc(sizeof(double) * particle_num_ * fitness_count);

	lower_bound_CPU = (double*)malloc(sizeof(double) * dim_);
	upper_bound_CPU = (double*)malloc(sizeof(double) * dim_);
	range_interval_CPU = (double*)malloc(sizeof(double) * dim_);
	//��ֵ
	for (int i = 0; i < dim_; i++)
	{
		lower_bound_CPU[i] = pso_para->lower_bound_[i];
		upper_bound_CPU[i] = pso_para->upper_bound_[i];
		range_interval_CPU[i] = upper_bound_CPU[i] - lower_bound_CPU[i];
	}


	//���Ӳ���GPU
	cudaMalloc((void**)& position_GPU, sizeof(double) * particle_num_ * dim_);
	cudaMalloc((void**)& velocity_GPU, sizeof(double) * particle_num_ * dim_);
	cudaMalloc((void**)& best_position_GPU, sizeof(double) * particle_num_ * dim_);
	cudaMalloc((void**)& fitness_GPU, sizeof(double) * particle_num_ * fitness_count);
	cudaMalloc((void**)& best_fitness_GPU, sizeof(double) * particle_num_ * fitness_count);


	cudaMalloc((void**)& upper_bound_, sizeof(double) * dim_);
	cudaMalloc((void**)& lower_bound_, sizeof(double) * dim_);
	cudaMalloc((void**)& range_interval_, sizeof(double) * dim_);
	//����CPU->GPU
	cudaMemcpy(upper_bound_, pso_para->upper_bound_, sizeof(double) * dim_, cudaMemcpyHostToDevice);
	cudaMemcpy(lower_bound_, pso_para->lower_bound_, sizeof(double) * dim_, cudaMemcpyHostToDevice);
	cudaMemcpy(range_interval_, pso_para->range_interval_, sizeof(double) * dim_, cudaMemcpyHostToDevice);


	cudaMalloc((void**)& all_best_position_, sizeof(double) * particle_num_ * dim_);
	cudaMalloc((void**)& all_best_fitness_, sizeof(double) * particle_num_ * fitness_count);

	//��ʼ�����������
	cudaMalloc(&globalState, particle_num_ * sizeof(curandState));
	initRandomGenerator << < 1, particle_num_ >> > (globalState, unsigned(time(NULL)));
	//��������������GPU�ռ�
	cudaMalloc((void**)& randomNumList, sizeof(double) * particle_num_);

	cudaMalloc((void**)& bestParticleIndex, sizeof(int));

	dt_ = pso_para->dt_;
	wstart_ = pso_para->wstart_;
	wend_ = pso_para->wend_;
	C1_ = pso_para->C1_;
	C2_ = pso_para->C2_;

	this->problemParas = problemParas;//CPU




	//��problemParas�Ĳ�����ֵ��GPU��
	DeviceSum = problemParas.DeviceSum;									//�豸����
	horiPointCount = problemParas.horiPointCount;						//δȥ��ǰ����ˮƽ����ĵ����Ŀ
	vertPointCount = problemParas.vertPointCount;						//δȥ��ǰ���д�ֱ����ĵ����Ŀ
	workShopLength = problemParas.workShopLength;						//���䳤��
	workShopWidth = problemParas.workShopWidth;							//������
	entrancePos = problemParas.entrancePos;								//�ֿ��������	
	exitPos = problemParas.exitPos;										//�ֿ��������

	//���ϲ����б�
	fixedLinkPointSum = 50;
	fixedUniqueLinkPointSum = 20;
	//���ͻ�����
	convey2DeviceDist = problemParas.convey2DeviceDist;//���ͻ����豸�ľ��루Ѱ·��ʱ��Ҫ���ǣ�
	conveyWidth = problemParas.conveyWidth;//���ͻ����
	conveyMinLength = problemParas.conveyMinLength;//���ͻ���̳���
	conveySpeed = problemParas.conveySpeed;//���ͻ������ٶ�
	strConveyorUnitCost = problemParas.strConveyorUnitCost;//��λֱ�����ͻ��ɱ�
	curveConveyorUnitCost = problemParas.curveConveyorUnitCost;//����ת�����ͻ��ɱ�
	conveyMinDist = problemParas.conveyMinDist;//������������֮�����̾���



	//DevicePara* deviceParaList;										//�豸�����б�
	//��Ŀ��DeviceSum
	cudaMalloc((void**)& ID, sizeof(int)* DeviceSum);
	cudaMalloc((void**)& workSpeed, sizeof(double)* DeviceSum);
	cudaMalloc((void**)& size, sizeof(Vector2)* DeviceSum);
	cudaMalloc((void**)& axis, sizeof(Vector2)* DeviceSum);
	cudaMalloc((void**)& direct, sizeof(DeviceDirect)* DeviceSum);
	cudaMalloc((void**)& spaceLength, sizeof(double)* DeviceSum);
	//����ڵ�����飨��Ӱ�������ߵĲ��֣�
	cudaMalloc((void**)& adjPInCount, sizeof(int)* DeviceSum);
	cudaMalloc((void**)& adjPOutCount, sizeof(int)* DeviceSum);
	cudaMalloc((void**)& accumAdjPInCount, sizeof(int)* DeviceSum);
	cudaMalloc((void**)& accumAdjPOutCount, sizeof(int)* DeviceSum);

	totalInPoint = problemParas.totalInPoint;
	totalOutPoint = problemParas.totalOutPoint;
	//cout << "totalInPoint:" << totalInPoint << endl;
	//cout << "totalOutPoint:" << totalOutPoint << endl;
	cudaMalloc((void**)& adjPointsIn, sizeof(AdjPoint)* totalInPoint);
	cudaMalloc((void**)& adjPointsOut, sizeof(AdjPoint)* totalOutPoint);



	//CPU��ֵ��GPU
	cudaMemcpy(ID, problemParas.ID, sizeof(int) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(workSpeed, problemParas.workSpeed, sizeof(double) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(size, problemParas.size, sizeof(Vector2) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(axis, problemParas.axis, sizeof(Vector2) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(direct, problemParas.direct, sizeof(DeviceDirect) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(spaceLength, problemParas.spaceLength, sizeof(double) * DeviceSum, cudaMemcpyHostToDevice);

	cudaMemcpy(adjPInCount, problemParas.adjPInCount, sizeof(int) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(adjPOutCount, problemParas.adjPOutCount, sizeof(int) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(accumAdjPInCount, problemParas.accumAdjPInCount, sizeof(int) * DeviceSum, cudaMemcpyHostToDevice);
	cudaMemcpy(accumAdjPOutCount, problemParas.accumAdjPOutCount, sizeof(int) * DeviceSum, cudaMemcpyHostToDevice);

	cudaMemcpy(adjPointsIn, problemParas.adjPointsIn, sizeof(int) * totalInPoint, cudaMemcpyHostToDevice);
	cudaMemcpy(adjPointsOut, problemParas.adjPointsOut, sizeof(int) * totalOutPoint, cudaMemcpyHostToDevice);


	//CargoType* cargoTypeList;						//���������б�
	//��Ŀ��CargoTypeNum
	CargoTypeNum = problemParas.CargoTypeNum;		//����������Ŀ
	totalLinkSum = problemParas.totalLinkSum;		//�ܵ���������Ŀ

	cudaMalloc((void**)& deviceSum, sizeof(int)* CargoTypeNum);
	cudaMalloc((void**)& linkSum, sizeof(int)* CargoTypeNum);
	cudaMalloc((void**)& accumLinkSum, sizeof(int)* CargoTypeNum);
	cudaMalloc((void**)& totalVolume, sizeof(double)* CargoTypeNum);
	cudaMalloc((void**)& deviceLinkList, sizeof(DeviceLink)* totalLinkSum);
	//CPU��ֵ��GPU
	cudaMemcpy(deviceSum, problemParas.deviceSum, sizeof(int)* CargoTypeNum, cudaMemcpyHostToDevice);
	cudaMemcpy(linkSum, problemParas.linkSum, sizeof(int)* CargoTypeNum, cudaMemcpyHostToDevice);
	cudaMemcpy(accumLinkSum, problemParas.accumLinkSum, sizeof(int)* CargoTypeNum, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceLinkList, problemParas.deviceLinkList, sizeof(DeviceLink)* totalLinkSum, cudaMemcpyHostToDevice);
	cudaMemcpy(totalVolume, problemParas.totalVolume, sizeof(double)* CargoTypeNum, cudaMemcpyHostToDevice);



	/*�������ӵ���������Ϣ�洢*/
	//��Ŀ��particle_num_
	//GPU�ڴ����
	cudaMalloc((void**)& curBestFitnessVal, sizeof(double) * particle_num_);
	inoutPSize = totalInPoint + totalOutPoint;
	//cudaMalloc((void**)& inoutPSize, sizeof(int) * particle_num_);
	cudaMalloc((void**)& inoutPoints, sizeof(InoutPoint) * inoutPSize * particle_num_);
	cudaMalloc((void**)& strConveyorList, sizeof(StraightConveyorInfo) * fixedUniqueLinkPointSum * totalLinkSum * particle_num_);
	cudaMalloc((void**)& strConveyorListSum, sizeof(int) * particle_num_);
	cudaMalloc((void**)& curveConveyorList, sizeof(Vector2Int) * fixedUniqueLinkPointSum * totalLinkSum * particle_num_);
	cudaMalloc((void**)& curveConveyorListSum, sizeof(int) * particle_num_);


	//��ǰ�������������Ϣ
	//��Ŀ1(�ȿ�����Ӧ��1��
	//GPU�ڴ����	
	cudaMalloc((void**)& curBestPath_FitnessVal, sizeof(double));
	cudaMalloc((void**)& curBestPath_InoutPoints, sizeof(InoutPoint) * inoutPSize);
	cudaMalloc((void**)& curBestPath_StrConveyorList, sizeof(StraightConveyorInfo)* fixedUniqueLinkPointSum * totalLinkSum);
	cudaMalloc((void**)& curBestPath_StrConveyorListSum, sizeof(int));
	cudaMalloc((void**)& curBestPath_CurveConveyorList, sizeof(Vector2Int) * fixedUniqueLinkPointSum * totalLinkSum);
	cudaMalloc((void**)& curBestPath_CurveConveyorListSum, sizeof(int));
	//��ʼ��
	InitCurBestPathFit << <1, 1 >> > (curBestPath_FitnessVal);


	cudaMalloc((void**)& pointDirectArray, sizeof(int) * 25);
	InitPointDirectArray << <1, 1 >> > (pointDirectArray);
}
//���µ�ǰ���������·��Ϣ
//ֻ��Ҫִ��һ��
static __global__ void UpdateCurBestPathInfo(int fitnessCount, double* fitness_GPU, int* bestParticleIndex,
	int tempStrConveyorList_PointSum, int tempCurveConveyorList_PointSum,
	/*���е���������Ϣ*/
	double* curBestFitnessVal, int inoutPSize, InoutPoint* inoutPoints, StraightConveyorInfo* strConveyorList,
	int* strConveyorListSum, Vector2Int* curveConveyorList, int* curveConveyorListSum,
	/*��ѵ���������Ϣ*/
	double* curBestPath_FitnessVal, int curBestPath_InoutPSize, InoutPoint* curBestPath_InoutPoints, StraightConveyorInfo* curBestPath_StrConveyorList,
	int* curBestPath_StrConveyorListSum, Vector2Int* curBestPath_CurveConveyorList, int* curBestPath_CurveConveyorListSum)
{
	//��curbest�����е���������ѱȽϣ�����curBest
	//���Ӧ��ֻҪ����һ��
	//ת������ֻ���ø��Ƶķ�����
	int bestIndex = bestParticleIndex[0];//������Ӷ�Ӧ���±�
	if (curBestPath_FitnessVal[0] > fitness_GPU[bestIndex * fitnessCount + 0]) {
		curBestPath_FitnessVal[0] = fitness_GPU[bestIndex * fitnessCount + 0];//ע��ƫ��ֵ
		for (int i = 0; i < curBestPath_InoutPSize; i++)
		{
			curBestPath_InoutPoints[i] = inoutPoints[bestIndex * inoutPSize + i];
		}
		for (int i = 0; i < strConveyorListSum[bestIndex]; i++) {
			curBestPath_StrConveyorList[i] = strConveyorList[bestIndex * tempStrConveyorList_PointSum + i];
		}
		curBestPath_StrConveyorListSum[0] = strConveyorListSum[bestIndex];
		for (int i = 0; i < curveConveyorListSum[bestIndex]; i++) {
			curBestPath_CurveConveyorList[i] = curveConveyorList[bestIndex * tempCurveConveyorList_PointSum + i];
		}
		curBestPath_CurveConveyorListSum[0] = curveConveyorListSum[bestIndex];
	}
}

PSOOptimizer::~PSOOptimizer()//������������Ҫ�޸�
{
	//CPU
	if (fitness_CPU) { delete[] fitness_CPU; }
	if (position_CPU) { delete[] position_CPU; }
	if (velocity_CPU) { delete[] velocity_CPU; }
	if (best_position_CPU) { delete[] best_position_CPU; }
	if (best_fitness_CPU) { delete[] best_fitness_CPU; }

	if (lower_bound_CPU) { delete[] lower_bound_CPU; }
	if (upper_bound_CPU) { delete[] upper_bound_CPU; }
	if (all_best_position_) { delete[] all_best_position_; }
	//GPU
	cudaFree(fitness_GPU);
	cudaFree(position_GPU);
	cudaFree(velocity_GPU);
	cudaFree(best_position_GPU);
	cudaFree(best_fitness_GPU);
	cudaFree(randomNumList);
	cudaFree(globalState);
	cudaFree(upper_bound_);
	cudaFree(lower_bound_);
	cudaFree(range_interval_);
	cudaFree(all_best_fitness_);
	cudaFree(all_best_position_);
	cudaFree(bestParticleIndex);

	//DevicePara* deviceParaList;
	cudaFree(ID);
	cudaFree(workSpeed);
	cudaFree(size);
	cudaFree(axis);
	cudaFree(direct);
	cudaFree(spaceLength);
	cudaFree(adjPInCount);
	cudaFree(adjPOutCount);
	cudaFree(accumAdjPInCount);
	cudaFree(accumAdjPOutCount);
	cudaFree(adjPointsIn);
	cudaFree(adjPointsOut);


	//CargoType* cargoTypeList;				//��������
	cudaFree(deviceSum);
	cudaFree(linkSum);
	cudaFree(accumLinkSum);
	cudaFree(deviceLinkList);
	cudaFree(totalVolume);


	//�������ӵ���������Ϣ
	cudaFree(curBestFitnessVal);
	cudaFree(inoutPoints);
	cudaFree(strConveyorList);
	cudaFree(strConveyorListSum);
	cudaFree(curveConveyorList);
	cudaFree(curveConveyorListSum);


	//��ǰ��ѵ���������Ϣ
	cudaFree(curBestPath_FitnessVal);
	cudaFree(curBestPath_InoutPoints);
	cudaFree(curBestPath_StrConveyorList);
	cudaFree(curBestPath_StrConveyorListSum);
	cudaFree(curBestPath_CurveConveyorList);
	cudaFree(curBestPath_CurveConveyorListSum);
	cudaFree(pointDirectArray);
}

// ��ʼ���������ӣ�û�и���ȫ����ѣ�
//CPU
void PSOOptimizer::InitialAllParticles()
{
	//��CPU�г�ʼ����������
	for (int i = 0; i < particle_num_; ++i)
	{
		InitialParticle(i);
	}

	//�������CPU�������ݿ���
	//for (int i = 0; i < particle_num_; i++)
	//{
	//	for (int j = 0; j < dim_; j++)
	//	{
	//		cout << position_CPU[i * dim_ + j] << ", ";
	//	}
	//	cout << endl;
	//}

	//���Ӵ�CPU->GPU
	cudaMemcpy(position_GPU, position_CPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyHostToDevice);
	cudaMemcpy(velocity_GPU, velocity_CPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyHostToDevice);
	cudaMemcpy(best_position_GPU, best_position_CPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyHostToDevice);
	cudaMemcpy(fitness_GPU, fitness_CPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyHostToDevice);
	cudaMemcpy(best_fitness_GPU, best_fitness_CPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyHostToDevice);
	//����Fitnessֵ GPU
	GetFitness();
	//��ʼ��pbest��gbest CPU
	for (int i = 0; i < particle_num_; ++i)
	{
		// ��ʼ����������λ��
		for (int j = 0; j < dim_; j++)
		{
			best_position_CPU[i * dim_ + j] = position_CPU[i * dim_ + j];
		}
		// ��ʼ�����Ӹ�����ʷ���
		for (int j = 0; j < fitness_count; j++)
		{
			best_fitness_CPU[i * fitness_count + j] = fitness_CPU[i * fitness_count + j];
		}
	}

}

// ��ʼ��Archive����
void PSOOptimizer::InitialArchiveList()
{
	//ͬ�����Ƚ���һ�����
	vector<Particle> particleList(particle_num_, Particle(dim_, fitness_count));
	for (int i = 0; i < particle_num_; i++)
	{
		//particleList[i].position_ = new double[dim_];
		//particleList[i].velocity_ = new double[dim_];
		//particleList[i].best_position_ = new double[dim_];
		for (int j = 0; j < dim_; j++)
		{
			particleList[i].position_[j] = position_CPU[i * dim_ + j];
			particleList[i].velocity_[j] = velocity_CPU[i * dim_ + j];
			particleList[i].best_position_[j] = best_position_CPU[i * dim_ + j];
		}

		//particleList[i].fitness_ = new double[fitness_count];
		//particleList[i].best_fitness_ = new double[fitness_count];
		for (int j = 0; j < fitness_count; j++)
		{
			particleList[i].fitness_[j] = fitness_CPU[i * fitness_count + j];
			particleList[i].best_fitness_[j] = best_fitness_CPU[i * fitness_count + j];
		}
	}
	//vector<Particle> particleList(this->particles_GPU, this->particles_GPU + this->particle_num_);
	Pareto initPareto(particleList);
	this->archive_list = initPareto.GetPareto();
}

// ����Archive���� �������ÿ�����Ӷ����㣬Ҳ����˵��ֻ��Ҫһ���߳̾Ϳ���
// Ϊ�˷��㣬����ʹ��CPU����
void PSOOptimizer::UpdateArchiveList()
{
	//GPU->CPU
	cudaMemcpy(position_CPU, position_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocity_CPU, velocity_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(best_position_CPU, best_position_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(fitness_CPU, fitness_GPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyDeviceToHost);
	cudaMemcpy(best_fitness_CPU, best_fitness_GPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyDeviceToHost);

	//���ȣ����㵱ǰ����Ⱥ��pareto�߽磬���߽����Ӽ��뵽�浵archiving��
	//����һ�����
	vector<Particle> particleList(particle_num_, Particle(dim_, fitness_count));
	for (int i = 0; i < particle_num_; i++)
	{
		//particleList[i].position_ = new double[dim_];
		//particleList[i].velocity_ = new double[dim_];
		//particleList[i].best_position_ = new double[dim_];
		for (int j = 0; j < dim_; j++)
		{
			particleList[i].position_[j] = position_CPU[i * dim_ + j];
			particleList[i].velocity_[j] = velocity_CPU[i * dim_ + j];
			particleList[i].best_position_[j] = best_position_CPU[i * dim_ + j];
		}

		//particleList[i].fitness_ = new double[fitness_count];
		//particleList[i].best_fitness_ = new double[fitness_count];
		for (int j = 0; j < fitness_count; j++)
		{
			particleList[i].fitness_[j] = fitness_CPU[i * fitness_count + j];
			particleList[i].best_fitness_[j] = best_fitness_CPU[i * fitness_count + j];
		}
	}
	//vector<Particle> particleList(this->particles_CPU, this->particles_CPU + particle_num_);
	Pareto pareto1(particleList);
	vector<Particle> curParetos = pareto1.GetPareto();
	//��Σ��ڴ浵�и���֧���ϵ���еڶ���ɸѡ�����Ǳ߽�����ȥ��
	vector<Particle> newParetos;
	curParetos.insert(curParetos.end(), this->archive_list.begin(), this->archive_list.end());//�ϲ�cur��ԭArchive
	Pareto pareto2(curParetos);
	vector<Particle> curArchives = pareto2.GetPareto();
	//���µ�ǰarchive
	this->archive_list = curArchives;
}

// ��ʼ��ȫ������
void PSOOptimizer::InitGbest()
{

	GetGbest getG(this->archive_list, this->meshDivCount, this->lower_bound_CPU, this->upper_bound_CPU, this->dim_, this->particle_num_);
	Particle* gbestList = getG.getGbest();
	//����gbestList�Ĳ������ݵ�CPU��Ȼ��GPU
	double* bestFitnessList = (double*)malloc(sizeof(double) * fitness_count * particle_num_);
	double* bestPositionList = (double*)malloc(sizeof(double) * dim_ * particle_num_);
	for (int i = 0; i < particle_num_; i++)
	{
		for (int j = 0; j < fitness_count; j++)
		{
			bestFitnessList[i * fitness_count + j] = gbestList[i].best_fitness_[j];
		}
		for (int k = 0; k < dim_; k++)
		{
			bestPositionList[i * dim_ + k] = gbestList[i].best_position_[k];
		}
	}
	//CPU->GPU
	cudaMemcpy(this->all_best_fitness_, bestFitnessList, sizeof(double) * fitness_count * particle_num_, cudaMemcpyHostToDevice);
	cudaMemcpy(this->all_best_position_, bestPositionList, sizeof(double) * dim_ * particle_num_, cudaMemcpyHostToDevice);
}

void PSOOptimizer::GetFitness()
{
	cout << "��ǰ������" << curr_iter_ << endl;
	FitnessFunction << <blockNum, threadsPerBlock >> > (curr_iter_, max_iter_num_, particle_num_, bestParticleIndex,
		/*ProblemParas proParas, �̶������ģ����ù�*/
		DeviceSum, fixedLinkPointSum, fixedUniqueLinkPointSum, vertPointCount, horiPointCount, workShopLength, workShopWidth, convey2DeviceDist, /*double conveyWidth, */
		strConveyorUnitCost, curveConveyorUnitCost, conveyMinDist, /*double conveyMinLength, */conveySpeed, entrancePos, exitPos,
		CargoTypeNum, totalLinkSum,

		/*CargoType* �̶�����*/
		linkSum, accumLinkSum, deviceLinkList, totalVolume,

		/*DevicePara**/
		size, spaceLength, adjPInCount, adjPOutCount, accumAdjPInCount, accumAdjPOutCount,
		totalInPoint, totalOutPoint, adjPointsIn, adjPointsOut,
		/*Particle*/
		dim_, fitness_count, fitness_GPU, position_GPU, /*double* velocity_GPU, double* best_positi on_GPU, double* best_fitness_GPU*/
		/*�洢��������������·��Ϣ*/
		curBestFitnessVal, inoutPSize, inoutPoints, strConveyorList, strConveyorListSum, curveConveyorList, curveConveyorListSum,
		pointDirectArray, globalState);
}

//�ȸĵ�һ��
void PSOOptimizer::UpdateAllParticles()
{
	//���㵱ǰ���Ĺ���ϵ��
	double temp = curr_iter_ / (double)max_iter_num_;
	//temp *= temp;//ϵ���仯
	w_ = wstart_ - (wstart_ - wend_) * temp;
	//���µ�ǰ����������
	//����ֻ�ܴ����洫��ȥ
	UpdateParticle_Kernal << <blockNum, threadsPerBlock >> > (curr_iter_, max_iter_num_, dim_, fitness_count, w_, C1_, C2_, dt_,
		fitness_GPU, position_GPU, velocity_GPU, best_position_GPU, best_fitness_GPU,
		globalState, randomNumList, range_interval_, upper_bound_, lower_bound_, all_best_position_,
		size, spaceLength, workShopLength, workShopWidth);//���allBest��CPU�ģ�
	GetFitness();
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		//goto Error;
	}
	curr_iter_++;
}

//����Pbest
void PSOOptimizer::UpdatePbest()
{
	//����GPU����
	UpdatePbest_Kernal << <blockNum, threadsPerBlock >> > (dim_, fitness_count, fitness_GPU, position_GPU, velocity_GPU,
		best_position_GPU, best_fitness_GPU, globalState);
}
// ����Gbest GPU
void PSOOptimizer::UpdateGbest()
{
	vector<Particle> tempArchiveL(this->archive_list);
	//GPU->CPU
	//cudaMemcpy(particles_CPU, particles_GPU, sizeof(Particle) * particle_num_, cudaMemcpyDeviceToHost);
	cudaMemcpy(position_CPU, position_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocity_CPU, velocity_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(best_position_CPU, best_position_GPU, sizeof(double) * particle_num_ * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(fitness_CPU, fitness_GPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyDeviceToHost);
	cudaMemcpy(best_fitness_CPU, best_fitness_GPU, sizeof(double) * particle_num_ * fitness_count, cudaMemcpyDeviceToHost);

	cudaMemcpy(lower_bound_CPU, lower_bound_, sizeof(double) * dim_, cudaMemcpyDeviceToHost);
	cudaMemcpy(upper_bound_CPU, upper_bound_, sizeof(double) * dim_, cudaMemcpyDeviceToHost);
	GetGbest getG(tempArchiveL, this->meshDivCount, lower_bound_CPU, upper_bound_CPU, this->dim_, this->particle_num_);
	Particle* gbestList = getG.getGbest();

	//GPU
	UpdateGbest_Kernal << <blockNum, threadsPerBlock >> > (fitness_count, dim_, all_best_fitness_, all_best_position_, gbestList);
}

//��ʼ����i������
//�����problemParas����ֱ����CPU��
void PSOOptimizer::InitialParticle(int i)
{
#pragma region ��ʼ��position/veloctiyֵ
	//���������Ȼ����ݳ���������ӵķ�Χ
	for (int j = 2; j < dim_; j += 3)
	{
		//cout << range_interval_CPU[j] << ", ";
		position_CPU[i * dim_ + j] = GetDoubleRand() * range_interval_CPU[j] + lower_bound_CPU[j];
		velocity_CPU[i * dim_ + j] = GetDoubleRand() * range_interval_CPU[j] / 300;
	}
	//cout << endl;
	//���ݳ����޸��豸���½緶Χ&�豸����
	vector<Vector2> deviceSizeCopy;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		deviceSizeCopy.push_back(Vector2(problemParas.size[i].x, problemParas.size[i].y));
	}
	for (int j = 2; j < dim_; j += 3)
	{
		//doubleתint��ת��ΪDirection��Ȼ����ݳ������¼����豸�ߴ�ͳ����
		//Rotate90����Rotate270,�޸�������
		DeviceDirect curDirect = (DeviceDirect)(int)position_CPU[i * dim_ + j];
		if (curDirect == DeviceDirect::Rotate90 || curDirect == DeviceDirect::Rotate270)
		{
			//x��y
			lower_bound_CPU[j - 2] = 0 + problemParas.size[j / 3].y * 0.5 + problemParas.spaceLength[j / 3];
			lower_bound_CPU[j - 1] = 0 + problemParas.size[j / 3].x * 0.5 + problemParas.spaceLength[j / 3];

			upper_bound_CPU[j - 2] = problemParas.workShopLength - problemParas.size[j / 3].y * 0.5 - problemParas.spaceLength[j / 3];
			upper_bound_CPU[j - 1] = problemParas.workShopWidth - problemParas.size[j / 3].x * 0.5 - problemParas.spaceLength[j / 3];

			//size��x��y��Ҫ����
			swap(deviceSizeCopy[j / 3].x, deviceSizeCopy[j / 3].y);
		}
		else
		{
			//x��y
			lower_bound_CPU[j - 2] = 0 + problemParas.size[j / 3].x * 0.5 + problemParas.spaceLength[j / 3];
			lower_bound_CPU[j - 1] = 0 + problemParas.size[j / 3].y * 0.5 + problemParas.spaceLength[j / 3];

			upper_bound_CPU[j - 2] = problemParas.workShopLength - problemParas.size[j / 3].x * 0.5 - problemParas.spaceLength[j / 3];
			upper_bound_CPU[j - 1] = problemParas.workShopWidth - problemParas.size[j / 3].y * 0.5 - problemParas.spaceLength[j / 3];
		}
		range_interval_CPU[j - 2] = upper_bound_CPU[j - 2] - lower_bound_CPU[j - 2];
		range_interval_CPU[j - 1] = upper_bound_CPU[j - 1] - lower_bound_CPU[j - 1];
	}


#pragma region ��ȫ���
	//for (int j = 0; j < dim_; j += 3) {
	//	particles_[i].position_[j] = GetDoubleRand() * range_interval_[j] + lower_bound_[j];
	//	particles_[i].position_[j + 1] = GetDoubleRand() * range_interval_[j + 1] + lower_bound_[j + 1];
	//	particles_[i].velocity_[j] = GetDoubleRand() * range_interval_[j] / 300;
	//	particles_[i].velocity_[j + 1] = GetDoubleRand() * range_interval_[j + 1] / 300;
	//}
#pragma endregion

#pragma region ���Ƿ��ص�Լ��������ֿ���������
//(ÿ��1�ײ���һ������㣬ֻҪ�ҵ�һ�������������ص�Լ�����Ͳ��ã�����Ĭ��Ϊ0
//�µ����������豸�İڷ�˳��
	vector<int> unmakeDeviceIndexVec;
	vector<int> madeDeviceIndexVec;
	for (int i = 0; i < problemParas.DeviceSum; i++)
	{
		unmakeDeviceIndexVec.push_back(i);
	}
	default_random_engine e;

	clock_t startTime, endTime;
	startTime = clock();//��ʱ��ʼ
	while (unmakeDeviceIndexVec.size() > 0)
	{
		endTime = clock();
		if ((static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC) > 1) {
#pragma region ��ȫ���
			for (int j = 0; j < dim_; j += 3) {
				position_CPU[i * dim_ + j] = GetDoubleRand() * range_interval_CPU[j] + lower_bound_CPU[j];
				position_CPU[i * dim_ + j + 1] = GetDoubleRand() * range_interval_CPU[j + 1] + lower_bound_CPU[j + 1];
				velocity_CPU[i * dim_ + j] = GetDoubleRand() * range_interval_CPU[j] / 300;
				velocity_CPU[i * dim_ + j + 1] = GetDoubleRand() * range_interval_CPU[j + 1] / 300;
			}
#pragma endregion
			break;
		}
		//΢�뼶���ȵ����������
		e.seed(GetRamdonSeed());
		uniform_int_distribution<unsigned> u(0, unmakeDeviceIndexVec.size() - 1);
		int randomVecIndex = u(e);
		int randomDeviceIndex = unmakeDeviceIndexVec[randomVecIndex];//�õ��豸��index
		int j = randomDeviceIndex * 3;

		double Xstart = lower_bound_CPU[j];
		double Ystart = lower_bound_CPU[j + 1];

		double tempPositionX = 0;
		double tempPositionY = 0;

		bool findParticle = false;
		while (Ystart <= upper_bound_CPU[j + 1] - 1 && findParticle == false) {//X��YҪ�ڷ�Χ��
			Xstart = lower_bound_CPU[j];
			while (Xstart <= upper_bound_CPU[j] - 1 && findParticle == false) {
				tempPositionX = GetDoubleRand() * 1.0 + Xstart;//�õ�Xstart��Xstart+1֮���һ�������
				tempPositionY = GetDoubleRand() * 1.0 + Ystart;//�õ�Ystart��Ystart+1֮���һ�������
				double halfX = deviceSizeCopy[j / 3].x * 0.5 + problemParas.spaceLength[j / 3];
				double halfY = deviceSizeCopy[j / 3].y * 0.5 + problemParas.spaceLength[j / 3];
				double tempLowX = tempPositionX - halfX;
				double tempUpX = tempPositionX + halfX;
				double tempLowY = tempPositionY - halfY;
				double tempUpY = tempPositionY + halfY;

				bool IsCross = false;
				//��鵱ǰ�豸�Ƿ��������ص�
				for (int k = 0; k < madeDeviceIndexVec.size(); k++)
				{
					int curDeviceIndex = madeDeviceIndexVec[k];
					int curDimIndex = curDeviceIndex * 3;

					double halfX1 = deviceSizeCopy[curDeviceIndex].x * 0.5 + problemParas.spaceLength[curDeviceIndex];
					double halfY1 = deviceSizeCopy[curDeviceIndex].y * 0.5 + problemParas.spaceLength[curDeviceIndex];

					double curLowX = position_CPU[i * dim_ + curDimIndex] - halfX1;
					double curUpX = position_CPU[i * dim_ + curDimIndex] + halfX1;
					double curLowY = position_CPU[i * dim_ + curDimIndex + 1] - halfY1;
					double curUpY = position_CPU[i * dim_ + curDimIndex + 1] + halfY1;
					//�������ص����˳�
					if (IsRangeOverlap(tempLowX, tempUpX, curLowX, curUpX) && IsRangeOverlap(tempLowY, tempUpY, curLowY, curUpY)) {
						IsCross = true;
						break;
					}
				}
				//ȫ�����ص��������Ӹ�ֵ
				if (IsCross == false) {
					findParticle = true;
					position_CPU[i * dim_ + j] = tempPositionX;
					position_CPU[i * dim_ + j + 1] = tempPositionY;
					velocity_CPU[i * dim_ + j] = GetDoubleRand() * range_interval_CPU[j] / 300;
					velocity_CPU[i * dim_ + j + 1] = GetDoubleRand() * range_interval_CPU[j + 1] / 300;

					//����vec
					madeDeviceIndexVec.push_back(randomDeviceIndex);
					unmakeDeviceIndexVec.erase(unmakeDeviceIndexVec.begin() + randomVecIndex);

				}
				Xstart++;
				if (Xstart >= upper_bound_CPU[j] - 1) {
					Ystart++;
				}
			}
		}
	}
#pragma endregion

#pragma endregion


}

