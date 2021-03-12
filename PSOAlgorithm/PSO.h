#pragma once
#include <vector>
#include <set>
#include <map>
#include <random>
#include "ProblemParas.h"

#include <cuda_runtime.h>//����cuda��API
#include <cuda_runtime_api.h>
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <curand.h>
#include <curand_kernel.h>
// ��Ӧ����Խ��Խ�û���ԽСԽ��
#define MINIMIZE_FITNESS
//#define MAXIMIZE_FITNESS

struct PSOPara
{
	int dim_;									// ����ά�ȣ�position��velocity��ά�ȣ�
	int fitness_count_;							// ��Ӧ����Ŀ
	int particle_num_;							// ���Ӹ���
	int max_iter_num_;							// ����������

	int mesh_div_count;

	double dt_;									// ʱ�䲽��
	double wstart_;								// ��ʼȨ��
	double wend_;								// ��ֹȨ��
	double C1_;									// ���ٶ�����1
	double C2_;									// ���ٶ�����2

	double* lower_bound_;						// position������Χ����
	double* upper_bound_;						// position������Χ����
	double* range_interval_;					// position�������䳤��

	int archive_max_count;						// pareto���Ž���������ֵ
	ProblemParas problemParas;					// �����Ӷ�Ӧ���豸���ֲ��� ȫ�ֲ���

	int blockSum;								// ÿ��Grid��block����Ŀ
	int threadsPerBlock;						// ÿ��Block��thread����Ŀ
	PSOPara() {}

	PSOPara(int dim)
	{
		dim_ = dim;

		lower_bound_ = new double[dim_];
		upper_bound_ = new double[dim_];
		range_interval_ = new double[dim_];

		ProblemParas problemParas();
	}

	// �����������ͷŶ��ڴ�
	~PSOPara()
	{
		if (lower_bound_) { delete[]lower_bound_; }
		if (upper_bound_) { delete[]upper_bound_; }
		if (range_interval_) { delete[]range_interval_; }
	}

	// ��������������ز���
	void SetProblemParas(ProblemParas paras) {
		problemParas = paras;
	}

	// ����ʱ�䲽��
	void SetDt(double stepLength) {
		dt_ = stepLength;
	}

	// ����wstart_
	void SetWstart(double startWeight) {
		wstart_ = startWeight;
	}

	// ����wend_
	void SetWend(double endWeight) {
		wend_ = endWeight;
	}

	// ����C1
	void SetC1(double c1) {
		C1_ = c1;
	}

	// ����C2
	void SetC2(double c2) {
		C2_ = c2;
	}

	// ����low_bound
	void SetLowBound(double lowBoundX, double lowBoundY, int lowBoundDirect) {
		for (int i = 0; i < dim_; i++) {
			if (i % 3 == 0)
				lower_bound_[i] = lowBoundX + problemParas.size[i / 3].x * 0.5 + problemParas.spaceLength[i / 3];
			else if (i % 3 == 1)
				lower_bound_[i] = lowBoundY + problemParas.size[i / 3].y * 0.5 + problemParas.spaceLength[i / 3];
			else
				lower_bound_[i] = lowBoundDirect;
		}
	}

	// ����upper_bound
	void SetUpBound(double upBoundX, double upBoundY, int upBoundDirect) {
		for (int i = 0; i < dim_; i++) {
			if (i % 3 == 0)
				upper_bound_[i] = upBoundX - problemParas.size[i / 3].x * 0.5 - problemParas.spaceLength[i / 3];
			else if (i % 3 == 1)
				upper_bound_[i] = upBoundY - problemParas.size[i / 3].y * 0.5 - problemParas.spaceLength[i / 3];
			else
				upper_bound_[i] = upBoundDirect;
		}
	}
};
//�洢�������ӵ���������Ϣ
struct BestPathInfo
{
	double curBestFitnessVal;//��ǰĿ�������ֵ�������ж��Ƿ���� fitnessCount
	int inoutPSize;//����ڵ������Ŀ	fitnessCount * 
	InoutPoint* inoutPoints;//����ڵ�ļ���  ��Ϊ���� 
	StraightConveyorInfo* strConveyorList;//ֱ�����ͻ���Ϣ�б� ��Ϊ����
	int strConveyorListSum;//fitnessCount
	Vector2Int* curveConveyorList;//ת�����ͻ���Ϣ�б� ��Ϊ����
	int curveConveyorListSum;//fitnessCount
	BestPathInfo() {
		curBestFitnessVal = INT_MAX;
	}
	void Clear() {

	}
};
//���ӽṹ��
struct Particle
{
	int dim_;							// ����ά�ȣ�position��velocity��ά�ȣ�
	int fitnessCount;					//��Ӧ�ȵĸ���
	double* fitness_;					//��Ӧ������
	double* position_;					//����λ������
	double* velocity_;					//�����ٶ�����

	double* best_position_;				//���ӵĸ�������λ������
	double* best_fitness_;				//���ӵĸ���������Ӧ������

	//vector<PointLink> pointLinks; //����·������

	//map<Vector2Int, PointInfo> pathPointInfoMap;//·������Ϣmap
	//set<SegPath> segPathSet;
	//int pointLinkSum = 0;//·������Ŀ
	Particle() {}
	Particle(int dim_, int fitnessCount)
	{
		this->dim_ = dim_;
		this->fitnessCount = fitnessCount;
		position_ = new double[dim_];
		velocity_ = new double[dim_];
		best_position_ = new double[dim_];
		fitness_ = new double[fitnessCount];
		best_fitness_ = new double[fitnessCount];
	}
	//���Ӧ�������ڱ����м�����ʱ��
	Particle(const Particle& particle)//�������캯��//Ҳ����˵��Ҫ���������ǳ����
	{
		this->dim_ = particle.dim_;
		this->fitnessCount = particle.fitnessCount;
		this->position_ = new double[this->dim_];
		this->velocity_ = new double[this->dim_];
		this->best_position_ = new double[this->dim_];

		this->fitness_ = new double[this->fitnessCount];
		this->best_fitness_ = new double[this->fitnessCount];
		for (int i = 0; i < this->fitnessCount; i++)
		{
			this->fitness_[i] = particle.fitness_[i];
			this->best_fitness_[i] = particle.best_fitness_[i];
		}
		for (int i = 0; i < this->dim_; i++)
		{
			this->position_[i] = particle.position_[i];
			this->velocity_[i] = particle.velocity_[i];
			this->best_position_[i] = particle.best_position_[i];
		}
	}
	~Particle()
	{

	}
};

class PSOOptimizer
{
public:
	int blockNum;							//block������Ŀ
	int threadsPerBlock;					//ÿ��block��thread��Ŀ

	int particle_num_;						//���Ӹ���	
	int max_iter_num_;						//����������
	int curr_iter_;							//��ǰ��������
	int dim_;								//����ά�ȣ�position��velocity��ά�ȣ�
	int fitness_count;						//��Ӧ����Ŀ
	int meshDivCount;						//����ȷ����ӣ�Ĭ��Ϊ10��
	int archiveMaxCount;					//archive����������Ŀ
	//Particle* particles_GPU;				//�������ӣ�GPU��
	//Particle* particles_CPU;				//�������ӣ�CPU��

	//���Ӳ���CPU
	double* fitness_CPU;					//��Ӧ������ �ܳ��ȣ�������Ŀ*fitnessCount
	double* position_CPU;					//����λ������ �ܳ��ȣ�������Ŀ*dim
	double* velocity_CPU;					//�����ٶ�����
	double* best_position_CPU;				//���ӵĸ�������λ������
	double* best_fitness_CPU;				//���ӵĸ���������Ӧ������

	double* lower_bound_CPU;				//position������Χ���� CPU
	double* upper_bound_CPU;				//position������Χ���� CPU
	double* range_interval_CPU;				//position�������䳤�� CPU
	//���Ӳ���GPU
	double* fitness_GPU;					//��Ӧ������
	double* position_GPU;					//����λ������
	double* velocity_GPU;					//�����ٶ�����
	double* best_position_GPU;				//���ӵĸ�������λ������
	double* best_fitness_GPU;				//���ӵĸ���������Ӧ������

	double* upper_bound_;					//position������Χ���� GPU
	double* lower_bound_;					//position������Χ���� GPU 
	double* range_interval_;				//position�������䳤�� GPU

	double* all_best_fitness_;				//ȫ���������ӵ���Ӧ������ 100x2 GPU
	double* all_best_position_;				//ȫ���������ӵ�position 100x12 GPU

	double* randomNumList;					//�������������
	curandState* globalState;				//cuda�����״̬����


	int* bestParticleIndex;					//ÿ�ֵ���������Fitness���ڵ�index

	double dt_;								//ʱ�䲽��
	double wstart_;							//��ʼȨ��
	double wend_;							//��ֹȨ��
	double w_;								//��ǰ����Ȩ��
	double C1_;								//���ٶ�����
	double C2_;								//���ٶ�����

	ProblemParas problemParas;				//����������� CPU


	//�����problemParas,GPU
	int DeviceSum;							//�豸����
	int horiPointCount;						//δȥ��ǰ����ˮƽ����ĵ����Ŀ
	int vertPointCount;						//δȥ��ǰ���д�ֱ����ĵ����Ŀ
	double workShopLength;					//���䳤��
	double workShopWidth;					//������
	Vector2 entrancePos;					//�ֿ��������	
	Vector2 exitPos;						//�ֿ��������

	//���ϲ����б�
	int fixedLinkPointSum;				//ÿһ��link�Ĺ̶�����ĿΪ50(û��ȥ��֮ǰ��)
	int fixedUniqueLinkPointSum;		//ȥ�غ��ÿһ��link�Ĺ̶�����ĿΪ20
	//���ͻ�����
	double convey2DeviceDist;				//���ͻ����豸�ľ��루Ѱ·��ʱ��Ҫ���ǣ�
	double conveyWidth;						//���ͻ����
	double conveyMinLength;					//���ͻ���̳���
	double conveySpeed;						//���ͻ������ٶ�
	double strConveyorUnitCost;				//��λֱ�����ͻ��ɱ�
	double curveConveyorUnitCost;			//����ת�����ͻ��ɱ�
	double conveyMinDist;					//������������֮�����̾���


	//DevicePara* deviceParaList;			//�豸�����б�
	//��Ŀ��DeviceSum
	int* ID;								//�豸ID
	double* workSpeed;						//�ӹ�/����1��λ���ϵ�ʱ��
	Vector2* size;							//�豸�ߴ磨�ֱ���x���y��ĳ��ȣ�
	Vector2* axis;							//�豸����
	DeviceDirect* direct;					//�豸����
	double* spaceLength;					//��϶��Ϊ��ʵ�־���Լ����
	//����ڵ�����飨��Ӱ�������ߵĲ��֣�
	int* adjPInCount;						//adjPointIn�������Ŀ����
	int* adjPOutCount;						//adjPointOut�������Ŀ����
	int* accumAdjPInCount;					//adjPointIn��Ŀ���ۼ�����
	int* accumAdjPOutCount;					//adjPointOut��Ŀ���ۼ�����

	int totalInPoint;						//��ڵ������Ŀ
	int totalOutPoint;						//���ڵ������Ŀ
	AdjPoint* adjPointsIn;					//��ڵ������
	AdjPoint* adjPointsOut;					//���ڵ������


	//CargoType* cargoTypeList;				//���������б�
	//��Ŀ��CargoTypeNum
	int CargoTypeNum;						//����������Ŀ
	int totalLinkSum;						//�ܵ���������Ŀ
	int* deviceSum;							//�������豸��Ŀ
	int* linkSum;							//�豸��Ե���Ŀ
	int* accumLinkSum;						//linkSum���ۼ�����
	DeviceLink* deviceLinkList;				//�豸�����б�
	double* totalVolume;					//�����ϵ���������


	vector<Particle> archive_list;			//���pareto���ӽ������ CPU

	//�洢ÿ�ֵ������������ӵ���������Ϣ
	//��ԭ���Ļ����� * ������Ŀ
	double* curBestFitnessVal;//��ǰĿ�������ֵ�������ж��Ƿ���� 1
	int inoutPSize;//����ڵ������Ŀ(totalInPoint + totalOutPoint)	1&��Ŀ�̶�
	InoutPoint* inoutPoints;//����ڵ�ļ���  totalInPoint + totalOutPoint

	StraightConveyorInfo* strConveyorList;//ֱ�����ͻ���Ϣ�б� fixedUniqueLinkPointSum * totalLinkSum
	int* strConveyorListSum;//ֱ�����ͻ���Ϣ��Ŀ�б���Ŀ���̶���
	Vector2Int* curveConveyorList;//ת�����ͻ���Ϣ�б� fixedUniqueLinkPointSum * totalLinkSum
	int* curveConveyorListSum;//ת�����ͻ���Ϣ��Ŀ�б�

	//�洢��ǰ���ŵ����ӵ���������Ϣ
	//��Ĭ��ֻ����Ӧ��1
	//BestPathInfo* bestPathInfoList;		//����·����Ϣ
	double* curBestPath_FitnessVal;//��ǰĿ�������ֵ�������ж��Ƿ���� 1
	InoutPoint* curBestPath_InoutPoints;//����ڵ�ļ���  totalInPoint + totalOutPoint
	StraightConveyorInfo* curBestPath_StrConveyorList;//ֱ�����ͻ���Ϣ�б� fixedUniqueLinkPointSum * totalLinkSum
	int* curBestPath_StrConveyorListSum;//ֱ�����ͻ���Ϣ��Ŀ�б� 1
	Vector2Int* curBestPath_CurveConveyorList;//ת�����ͻ���Ϣ�б� fixedUniqueLinkPointSum * totalLinkSum
	int* curBestPath_CurveConveyorListSum;//ת�����ͻ���Ϣ��Ŀ�б� 1


	//Up = 1, Right = 2, Down = 3, Left = 4
	//һ��10����������£����ң����£����ϣ��������ң�
	//�������ң���������
	//������˳������������
	int* pointDirectArray;
public:
	// Ĭ�Ϲ��캯��
	PSOOptimizer() {}

	// ���캯��
	PSOOptimizer(PSOPara* pso_para, ProblemParas& problemParas);

	// ��������
	~PSOOptimizer();

	// ��ʼ���������Ӳ���
	void InitialAllParticles();

	// ��ʼ����i�����Ӳ���
	void InitialParticle(int i);

	// ��ʼ��Archive����
	void InitialArchiveList();

	// ����Archive����
	void UpdateArchiveList();

	// ��ʼ��ȫ������
	void InitGbest();

	// ��������ӵ���Ӧ��ֵ
	void GetFitness();

	// �����������Ӳ���
	void UpdateAllParticles();

	// ���µ�i������
	//void UpdateParticle(int i);

	// ����Pbest
	void UpdatePbest();

	// ����Gbest
	void UpdateGbest();

	// �Ƚ��������ӵ���Ӧ�ȣ��ж��Ƿ���ȫ֧�䣬�Ӷ������pbest
	bool ComparePbest(double* fitness, double* pbestFitness);


};