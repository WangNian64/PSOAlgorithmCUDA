#pragma once

#include "PSO.h"
class MeshCrowd
{
public:
	vector<Particle> curArchiveList;	//存档中的所有粒子
	int meshDivCount;			//网格等分因子，默认为10
	int curArchiveLength;		//存档中粒子数目
	int ParticleCount;			//粒子群数目

	int* archiveMeshIdList;		//粒子编号
	int* archiveCrowdList;		//拥挤度矩阵

	Particle* gbestParticleList;//gbest粒子数组

	int dim;
	double* minList;
	double* maxList;

public:
	// 默认构造函数
	MeshCrowd() {}

	// 构造函数
	MeshCrowd(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount);

	// 析构函数
	~MeshCrowd();

	// 计算网格id
	int CalMeshId(Particle particle);

	// 对每个粒子定义网格编号
	void DivideArchiving();

	// 计算拥挤度列表
	void GetCrowd();
};

class GetGbest : MeshCrowd
{
public:
	double* archiveProbability;

	// 默认构造函数
	GetGbest() {}

	// 构造函数
	GetGbest(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount);

	// 计算粒子被选择的概率
	void GetProbability();

	// 
	int GetGbestIndex();

	// 计算当前粒子群中每一个粒子的全局最优解集
	Particle* getGbest();



	// 计算粒子被选择的概率
	void GetProbability1();

	// 按概率清除粒子，拥挤度高的粒子被清除的概率越高
	vector<int> GetClearIndex(int threshN);

	// Clear方法
	vector<Particle> Clear(int thresh);
};

