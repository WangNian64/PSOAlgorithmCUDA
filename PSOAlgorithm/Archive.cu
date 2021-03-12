#pragma once
#include "Archive.h"
// ���캯��
MeshCrowd::MeshCrowd(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount)
{
	cout << "����MeshCrowd" << endl;
	this->curArchiveList = curArchiveList;//�浵��������
	this->curArchiveLength = curArchiveList.size();//�浵����

	this->meshDivCount = meshDivCount;//����ȷ����ӣ�Ĭ��Ϊ10
	this->ParticleCount = ParticleCount;//����Ⱥ����

	this->archiveMeshIdList = new int[this->curArchiveLength];//���ӱ��
	this->archiveCrowdList = new int[this->curArchiveLength];//ӵ���Ⱦ������ڼ�¼��ǰ�������������
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		this->archiveMeshIdList[i] = 0;
		this->archiveCrowdList[i] = 0;
	}
	this->gbestParticleList = new Particle[this->ParticleCount];//��ʼ��gbest�����б�

	this->minList = new double[dim];
	this->maxList = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		this->minList[i] = minList[i];
		this->maxList[i] = maxList[i];
	}
	cout << "����MeshCrowd" << endl;
}

// ��������
MeshCrowd::~MeshCrowd()
{
	if (archiveMeshIdList) { delete[]archiveMeshIdList; }
	if (archiveCrowdList) { delete[]archiveCrowdList; }

	if (gbestParticleList) { delete[]gbestParticleList; }
	if (minList) { delete[]minList; }
	if (maxList) { delete[]maxList; }
}

// ��������id
int MeshCrowd::CalMeshId(Particle particle)
{
	//����������id
	//���ȣ���ÿ��ά�Ȱ��յȷ����ӽ��еȷ���ɢ����
	//��ȡ�����ڸ�ά���ϵı�š�����10���ƽ�ÿһ��ά�ȱ�ŵȱ���ӣ�����û��Զ�����mesh_div_num��ֵ�������Զ��壩,�����ֵ
	//Ҳ���ǣ�ÿ��ά���϶���һ��ֵ
	int id = 0;
	for (int i = 0; i < particle.dim_; i++)
	{
		int id_dim = (int)((particle.position_[i] - minList[i]) * meshDivCount / (maxList[i] - minList[i]));
		id = id + id_dim * (pow(meshDivCount, i));
	}
	return id;
}

// ��ÿ�����Ӷ��������ţ������е����ӷŵ���ͬ�������У�
void MeshCrowd::DivideArchiving()
{
	for (int i = 0; i < this->curArchiveLength; i++)
	{
		this->archiveMeshIdList[i] = CalMeshId(curArchiveList[i]);
	}
}

// �����������������ӵ�ӵ����
void MeshCrowd::GetCrowd()
{
	//����һ�����������Ӽ��������ţ����ڸ�������
	vector<int> particleIndexList(this->curArchiveLength);
	//��ʼ��indexlist��ֵ
	for (int i = 0; i < particleIndexList.size(); i++)
	{
		particleIndexList[i] = i;
	}
	while (particleIndexList.size() > 0)
	{
		//��ű�����ѭ������index[0]���Ӿ�����ͬ����id���м���λ
		vector<int> index_same;
		index_same.push_back(particleIndexList[0]);
		for (int i = 1; i < particleIndexList.size(); i++)
		{
			if (this->archiveMeshIdList[particleIndexList[0]] == this->archiveMeshIdList[particleIndexList[i]])
			{
				index_same.push_back(particleIndexList[i]);
			}
		}
		int particleNum = index_same.size();//����������������
		for (int i = 0; i < index_same.size(); i++)//���±�������id�µ��������ӵ�ӵ����
		{
			this->archiveCrowdList[index_same[i]] = particleNum;
			//�ҵ�particleIndexList��index_same[i]��ƥ��λ��
			auto it = find(particleIndexList.begin(), particleIndexList.end(), index_same[i]);
			particleIndexList.erase(it);//ɾ���������������������Ӷ�Ӧ�������ţ������ظ�����
		}
	}

}

// ���캯��
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

// �������ӱ�ѡ��ĸ���
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
		archiveProbability[i] = archiveProbability[i] / totalProb;//��һ�������ӱ�ѡ������ۺ�Ϊ1
	}
}

// ��������
int GetGbest::GetGbestIndex()
{
	double randomProb = rand() % 1000 / (double)1000;//�������С��
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

// ���㵱ǰ����Ⱥ��ÿһ�����ӵ�ȫ�����Ž⼯
Particle* GetGbest::getGbest()
{
	this->GetProbability();//����浵��ÿ�����Ӻ�ӵ������صĸ���
	//����ӵ���ȸߵ����ѡ��gbest
	for (int i = 0; i < this->ParticleCount; i++)
	{
		int gbestIndex = this->GetGbestIndex();
		this->gbestParticleList[i] = this->curArchiveList[gbestIndex];
	}
	return this->gbestParticleList;
}


// �������ӱ�ѡ��ĸ���1
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

// ������������ӣ�ӵ���ȸߵ����ӱ�����ĸ���Խ��
vector<int> GetGbest::GetClearIndex(int thresh)
{
	int clearLen = this->curArchiveList.size() - thresh;//��Ҫ���������������
	vector<int> clearIndexList;
	while (clearIndexList.size() < clearLen)
	{
		double randomProb = rand() % 1000 / (double)1000;//�������С��
		for (int i = 0; i < this->curArchiveLength; i++)
		{
			double subTotal = 0.0;
			for (int j = 0; j < i + 1; j++)
			{
				subTotal += this->archiveProbability[j];
			}
			if (randomProb <= subTotal)
			{
				//����i�Ĳ��ӽ�ȥ
				if (find(clearIndexList.begin(), clearIndexList.end(), i) == clearIndexList.end())
				{
					clearIndexList.push_back(i);//��¼����ֵ
					break;
				}
			}
		}
	}
	return clearIndexList;
}

// Clear����
vector<Particle> GetGbest::Clear(int thresh)
{
	this->GetProbability1();
	vector<int> clearIndexList = this->GetClearIndex(thresh);//��Ҫɾ�����±�����
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