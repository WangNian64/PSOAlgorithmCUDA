#pragma once
#include <vector>
#include <set>
#include <map>
#include <random>
#include "ProblemParas.h"

#include <cuda_runtime.h>//几个cuda的API
#include <cuda_runtime_api.h>
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <curand.h>
#include <curand_kernel.h>
// 适应度是越大越好还是越小越好
#define MINIMIZE_FITNESS
//#define MAXIMIZE_FITNESS

struct PSOPara
{
	int dim_;									// 参数维度（position和velocity的维度）
	int fitness_count_;							// 适应度数目
	int particle_num_;							// 粒子个数
	int max_iter_num_;							// 最大迭代次数

	int mesh_div_count;

	double dt_;									// 时间步长
	double wstart_;								// 初始权重
	double wend_;								// 终止权重
	double C1_;									// 加速度因子1
	double C2_;									// 加速度因子2

	double* lower_bound_;						// position搜索范围下限
	double* upper_bound_;						// position搜索范围上限
	double* range_interval_;					// position搜索区间长度

	int archive_max_count;						// pareto最优解数组的最大值
	ProblemParas problemParas;					// 和粒子对应的设备布局参数 全局参数

	int blockSum;								// 每个Grid中block的数目
	int threadsPerBlock;						// 每个Block中thread的数目
	PSOPara() {}

	PSOPara(int dim)
	{
		dim_ = dim;

		lower_bound_ = new double[dim_];
		upper_bound_ = new double[dim_];
		range_interval_ = new double[dim_];

		ProblemParas problemParas();
	}

	// 析构函数：释放堆内存
	~PSOPara()
	{
		if (lower_bound_) { delete[]lower_bound_; }
		if (upper_bound_) { delete[]upper_bound_; }
		if (range_interval_) { delete[]range_interval_; }
	}

	// 设置物流问题相关参数
	void SetProblemParas(ProblemParas paras) {
		problemParas = paras;
	}

	// 设置时间步长
	void SetDt(double stepLength) {
		dt_ = stepLength;
	}

	// 设置wstart_
	void SetWstart(double startWeight) {
		wstart_ = startWeight;
	}

	// 设置wend_
	void SetWend(double endWeight) {
		wend_ = endWeight;
	}

	// 设置C1
	void SetC1(double c1) {
		C1_ = c1;
	}

	// 设置C2
	void SetC2(double c2) {
		C2_ = c2;
	}

	// 设置low_bound
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

	// 设置upper_bound
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
//存储最优粒子的输送线信息
struct BestPathInfo
{
	double curBestFitnessVal;//当前目标的最优值，用于判断是否更新 fitnessCount
	int inoutPSize;//出入口点的总数目	fitnessCount * 
	InoutPoint* inoutPoints;//出入口点的集合  人为设置 
	StraightConveyorInfo* strConveyorList;//直线输送机信息列表 人为设置
	int strConveyorListSum;//fitnessCount
	Vector2Int* curveConveyorList;//转弯输送机信息列表 人为设置
	int curveConveyorListSum;//fitnessCount
	BestPathInfo() {
		curBestFitnessVal = INT_MAX;
	}
	void Clear() {

	}
};
//粒子结构体
struct Particle
{
	int dim_;							// 参数维度（position和velocity的维度）
	int fitnessCount;					//适应度的个数
	double* fitness_;					//适应度数组
	double* position_;					//粒子位置数组
	double* velocity_;					//粒子速度数组

	double* best_position_;				//粒子的个体最优位置数组
	double* best_fitness_;				//粒子的个体最优适应度数组

	//vector<PointLink> pointLinks; //最终路径集合

	//map<Vector2Int, PointInfo> pathPointInfoMap;//路径点信息map
	//set<SegPath> segPathSet;
	//int pointLinkSum = 0;//路径的数目
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
	//深拷贝应该是用在保存中间结果的时候
	Particle(const Particle& particle)//拷贝构造函数//也就是说需要深拷贝而不是浅拷贝
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
	int blockNum;							//block的总数目
	int threadsPerBlock;					//每个block的thread数目

	int particle_num_;						//粒子个数	
	int max_iter_num_;						//最大迭代次数
	int curr_iter_;							//当前迭代次数
	int dim_;								//参数维度（position和velocity的维度）
	int fitness_count;						//适应度数目
	int meshDivCount;						//网格等分因子（默认为10）
	int archiveMaxCount;					//archive数组的最大数目
	//Particle* particles_GPU;				//所有粒子（GPU）
	//Particle* particles_CPU;				//所有粒子（CPU）

	//粒子参数CPU
	double* fitness_CPU;					//适应度数组 总长度：粒子数目*fitnessCount
	double* position_CPU;					//粒子位置数组 总长度：粒子数目*dim
	double* velocity_CPU;					//粒子速度数组
	double* best_position_CPU;				//粒子的个体最优位置数组
	double* best_fitness_CPU;				//粒子的个体最优适应度数组

	double* lower_bound_CPU;				//position搜索范围上限 CPU
	double* upper_bound_CPU;				//position搜索范围下限 CPU
	double* range_interval_CPU;				//position搜索区间长度 CPU
	//粒子参数GPU
	double* fitness_GPU;					//适应度数组
	double* position_GPU;					//粒子位置数组
	double* velocity_GPU;					//粒子速度数组
	double* best_position_GPU;				//粒子的个体最优位置数组
	double* best_fitness_GPU;				//粒子的个体最优适应度数组

	double* upper_bound_;					//position搜索范围上限 GPU
	double* lower_bound_;					//position搜索范围下限 GPU 
	double* range_interval_;				//position搜索区间长度 GPU

	double* all_best_fitness_;				//全局最优粒子的适应度数组 100x2 GPU
	double* all_best_position_;				//全局最优粒子的position 100x12 GPU

	double* randomNumList;					//存随机数的数组
	curandState* globalState;				//cuda随机数状态数组


	int* bestParticleIndex;					//每轮迭代中最优Fitness所在的index

	double dt_;								//时间步长
	double wstart_;							//初始权重
	double wend_;							//终止权重
	double w_;								//当前迭代权重
	double C1_;								//加速度因子
	double C2_;								//加速度因子

	ProblemParas problemParas;				//布局问题参数 CPU


	//拆解后的problemParas,GPU
	int DeviceSum;							//设备总数
	int horiPointCount;						//未去重前所有水平方向的点的数目
	int vertPointCount;						//未去重前所有垂直方向的点的数目
	double workShopLength;					//车间长度
	double workShopWidth;					//车间宽度
	Vector2 entrancePos;					//仓库入口坐标	
	Vector2 exitPos;						//仓库出口坐标

	//物料参数列表
	int fixedLinkPointSum;				//每一条link的固定点数目为50(没有去重之前的)
	int fixedUniqueLinkPointSum;		//去重后的每一条link的固定点数目为20
	//输送机参数
	double convey2DeviceDist;				//输送机到设备的距离（寻路的时候要考虑）
	double conveyWidth;						//输送机宽度
	double conveyMinLength;					//输送机最短长度
	double conveySpeed;						//输送机输送速度
	double strConveyorUnitCost;				//单位直线输送机成本
	double curveConveyorUnitCost;			//单个转弯输送机成本
	double conveyMinDist;					//输送线两个点之间的最短距离


	//DevicePara* deviceParaList;			//设备参数列表
	//数目：DeviceSum
	int* ID;								//设备ID
	double* workSpeed;						//加工/处理1单位物料的时间
	Vector2* size;							//设备尺寸（分别是x轴和y轴的长度）
	Vector2* axis;							//设备坐标
	DeviceDirect* direct;					//设备朝向
	double* spaceLength;					//空隙（为了实现距离约束）
	//出入口点的数组（会影响输送线的布局）
	int* adjPInCount;						//adjPointIn数组的数目数组
	int* adjPOutCount;						//adjPointOut数组的数目数组
	int* accumAdjPInCount;					//adjPointIn数目的累加数组
	int* accumAdjPOutCount;					//adjPointOut数目的累加数组

	int totalInPoint;						//入口点的总数目
	int totalOutPoint;						//出口点的总数目
	AdjPoint* adjPointsIn;					//入口点的数组
	AdjPoint* adjPointsOut;					//出口点的数组


	//CargoType* cargoTypeList;				//货物类型列表
	//数目：CargoTypeNum
	int CargoTypeNum;						//货物类型数目
	int totalLinkSum;						//总的连接线数目
	int* deviceSum;							//经过的设备数目
	int* linkSum;							//设备配对的数目
	int* accumLinkSum;						//linkSum的累加数组
	DeviceLink* deviceLinkList;				//设备连接列表
	double* totalVolume;					//该物料的总物流量


	vector<Particle> archive_list;			//存放pareto非劣解的数组 CPU

	//存储每轮迭代中所有粒子的输送线信息
	//在原来的基础上 * 粒子数目
	double* curBestFitnessVal;//当前目标的最优值，用于判断是否更新 1
	int inoutPSize;//出入口点的总数目(totalInPoint + totalOutPoint)	1&数目固定
	InoutPoint* inoutPoints;//出入口点的集合  totalInPoint + totalOutPoint

	StraightConveyorInfo* strConveyorList;//直线输送机信息列表 fixedUniqueLinkPointSum * totalLinkSum
	int* strConveyorListSum;//直线输送机信息数目列表（数目不固定）
	Vector2Int* curveConveyorList;//转弯输送机信息列表 fixedUniqueLinkPointSum * totalLinkSum
	int* curveConveyorListSum;//转弯输送机信息数目列表

	//存储当前最优的粒子的输送线信息
	//先默认只存适应度1
	//BestPathInfo* bestPathInfoList;		//最优路径信息
	double* curBestPath_FitnessVal;//当前目标的最优值，用于判断是否更新 1
	InoutPoint* curBestPath_InoutPoints;//出入口点的集合  totalInPoint + totalOutPoint
	StraightConveyorInfo* curBestPath_StrConveyorList;//直线输送机信息列表 fixedUniqueLinkPointSum * totalLinkSum
	int* curBestPath_StrConveyorListSum;//直线输送机信息数目列表 1
	Vector2Int* curBestPath_CurveConveyorList;//转弯输送机信息列表 fixedUniqueLinkPointSum * totalLinkSum
	int* curBestPath_CurveConveyorListSum;//转弯输送机信息数目列表 1


	//Up = 1, Right = 2, Down = 3, Left = 4
	//一共10种情况：上下，左右，下下，上上，左左，右右，
	//上左，上右，下左，下右
	//横竖的顺序都是上下左右
	int* pointDirectArray;
public:
	// 默认构造函数
	PSOOptimizer() {}

	// 构造函数
	PSOOptimizer(PSOPara* pso_para, ProblemParas& problemParas);

	// 析构函数
	~PSOOptimizer();

	// 初始化所有粒子参数
	void InitialAllParticles();

	// 初始化第i个粒子参数
	void InitialParticle(int i);

	// 初始化Archive数组
	void InitialArchiveList();

	// 更新Archive数组
	void UpdateArchiveList();

	// 初始化全局最优
	void InitGbest();

	// 计算该粒子的适应度值
	void GetFitness();

	// 更新所有粒子参数
	void UpdateAllParticles();

	// 更新第i个粒子
	//void UpdateParticle(int i);

	// 更新Pbest
	void UpdatePbest();

	// 更新Gbest
	void UpdateGbest();

	// 比较两个粒子的适应度，判断是否完全支配，从而计算出pbest
	bool ComparePbest(double* fitness, double* pbestFitness);


};