#pragma once
#include "Archive.h"
// 构造函数
MeshCrowd::MeshCrowd(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount)
{
	cout << "进入MeshCrowd" << endl;
	this->curArchiveList = curArchiveList;//存档粒子数组
	this->curArchiveLength = curArchiveList.size();//存档长度

	this->meshDivCount = meshDivCount;//网格等分因子，默认为10
	this->ParticleCount = ParticleCount;//粒子群数量

	this->archiveMeshIdList = new int[this->curArchiveLength];//粒子编号
	this->archiveCrowdList = new int[this->curArchiveLength];//拥挤度矩阵，用于记录当前网格的粒子数量
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		this->archiveMeshIdList[i] = 0;
		this->archiveCrowdList[i] = 0;
	}
	this->gbestParticleList = new Particle[this->ParticleCount];//初始化gbest粒子列表

	this->minList = new double[dim];
	this->maxList = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		this->minList[i] = minList[i];
		this->maxList[i] = maxList[i];
	}
	cout << "进入MeshCrowd" << endl;
}

// 析构函数
MeshCrowd::~MeshCrowd()
{
	if (archiveMeshIdList) { delete[]archiveMeshIdList; }
	if (archiveCrowdList) { delete[]archiveCrowdList; }

	if (gbestParticleList) { delete[]gbestParticleList; }
	if (minList) { delete[]minList; }
	if (maxList) { delete[]maxList; }
}

// 计算网格id
int MeshCrowd::CalMeshId(Particle particle)
{
	//计算网格编号id
	//首先，将每个维度按照等分因子进行等分离散化，
	//获取粒子在各维度上的编号。按照10进制将每一个维度编号等比相加（如果用户自定义了mesh_div_num的值，则按照自定义）,计算出值
	//也就是，每个维度上都有一个值
	int id = 0;
	for (int i = 0; i < particle.dim_; i++)
	{
		int id_dim = (int)((particle.position_[i] - minList[i]) * meshDivCount / (maxList[i] - minList[i]));
		id = id + id_dim * (pow(meshDivCount, i));
	}
	return id;
}

// 对每个粒子定义网格编号（将所有的粒子放到不同的网格中）
void MeshCrowd::DivideArchiving()
{
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		this->archiveMeshIdList[i] = CalMeshId(curArchiveList[i]);
	}
}

// 计算所有网格中粒子的拥挤度
void MeshCrowd::GetCrowd()
{
	//定义一个数组存放粒子集的索引号，用于辅助计算
	vector<int> particleIndexList(this->curArchiveLength);
	//初始化indexlist的值
	for (int i = 0; i < particleIndexList.size(); i++)
	{
		particleIndexList[i] = i;
	}
	while (particleIndexList.size() > 0)
	{
		//存放本次子循环中与index[0]粒子具有相同网格id所有检索位
		vector<int> index_same;
		index_same.push_back(particleIndexList[0]);
		for (int i = 1; i < particleIndexList.size(); i++)
		{
			if (this->archiveMeshIdList[particleIndexList[0]] == this->archiveMeshIdList[particleIndexList[i]])
			{
				index_same.push_back(particleIndexList[i]);
			}
		}
		int particleNum = index_same.size();//本轮网格中粒子数
		for (int i = 0; i < index_same.size(); i++)//更新本轮网格id下的所有粒子的拥挤度
		{
			this->archiveCrowdList[index_same[i]] = particleNum;
			//找到particleIndexList中index_same[i]的匹配位置
			auto it = find(particleIndexList.begin(), particleIndexList.end(), index_same[i]);
			particleIndexList.erase(it);//删除本轮网格所包含的粒子对应的索引号，避免重复计算
		}
	}

}

// 构造函数
GetGbest::GetGbest(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount) :
	MeshCrowd(curArchiveList, meshDivCount, minList, maxList, dim, ParticleCount)
{
	archiveProbability = new double[this->curArchiveLength];
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		archiveProbability[i] = 0;
	}
	this->DivideArchiving();
	this->GetCrowd();
}

// 计算粒子被选择的概率
void GetGbest::GetProbability()
{
	double totalProb = 0.0;
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		archiveProbability[i] = 1.0 / pow(archiveCrowdList[i], 3);
		totalProb += archiveProbability[i];
	}
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		archiveProbability[i] = archiveProbability[i] / totalProb;//归一化，粒子被选择概率综合为1
	}
}

// 返回索引
int GetGbest::GetGbestIndex()
{
	double randomProb = rand() % 1000 / (double)1000;//产生随机小数
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		double subProbTotal = 0.0;
		for (int j = 0; j < i + 1; j++)
		{
			subProbTotal += this->archiveProbability[j];
		}
		if (randomProb <= subProbTotal)
		{
			return i;
		}
	}
}

// 计算当前粒子群中每一个粒子的全局最优解集
Particle* GetGbest::getGbest()
{
	this->GetProbability();//计算存档中每个粒子和拥挤度相关的概率
	//按照拥挤度高低随机选择gbest
	for (int i = 0; i < this->ParticleCount; i++)
	{
		int gbestIndex = this->GetGbestIndex();
		this->gbestParticleList[i] = this->curArchiveList[gbestIndex];
	}
	return this->gbestParticleList;
}


// 计算粒子被选择的概率1
void GetGbest::GetProbability1()
{
	int total = 0;
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		total += this->archiveCrowdList[i];
	}
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		this->archiveProbability[i] = this->archiveCrowdList[i] / total;
	}
}

// 按概率清除粒子，拥挤度高的粒子被清除的概率越高
vector<int> GetGbest::GetClearIndex(int thresh)
{
	int clearLen = this->curArchiveList.size() - thresh;//需要清除掉的粒子数量
	vector<int> clearIndexList;
	while (clearIndexList.size() < clearLen)
	{
		double randomProb = rand() % 1000 / (double)1000;//产生随机小数
		for (int i = 0; i < this->curArchiveLength; i++)
		{
			double subTotal = 0.0;
			for (int j = 0; j < i + 1; j++)
			{
				subTotal += this->archiveProbability[j];
			}
			if (randomProb <= subTotal)
			{
				//等于i的不加进去
				if (find(clearIndexList.begin(), clearIndexList.end(), i) == clearIndexList.end())
				{
					clearIndexList.push_back(i);//记录检索值
					break;
				}
			}
		}
	}
	return clearIndexList;
}

// Clear方法
vector<Particle> GetGbest::Clear(int thresh)
{
	this->GetProbability1();
	vector<int> clearIndexList = this->GetClearIndex(thresh);//需要删除的下标数组
	vector<Particle> resultArchive;
	int curIndex = 0;
	for (auto iter = curArchiveList.begin(); iter != curArchiveList.end(); iter++)
	{
		if (curIndex < clearIndexList.size())
		{
			if (iter != curArchiveList.begin() + clearIndexList[curIndex])
			{
				resultArchive.push_back(*iter);
			}
			else
			{
				curIndex++;
			}
		}
	}
	return resultArchive;
}