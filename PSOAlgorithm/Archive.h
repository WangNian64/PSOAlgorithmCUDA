#pragma once

#include "PSO.h"
class MeshCrowd
{
public:
	vector<Particle> curArchiveList;	//�浵�е���������
	int meshDivCount;			//����ȷ����ӣ�Ĭ��Ϊ10
	int curArchiveLength;		//�浵��������Ŀ
	int ParticleCount;			//����Ⱥ��Ŀ

	int* archiveMeshIdList;		//���ӱ��
	int* archiveCrowdList;		//ӵ���Ⱦ���

	Particle* gbestParticleList;//gbest��������

	int dim;
	double* minList;
	double* maxList;

public:
	// Ĭ�Ϲ��캯��
	MeshCrowd() {}

	// ���캯��
	MeshCrowd(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount);

	// ��������
	~MeshCrowd();

	// ��������id
	int CalMeshId(Particle particle);

	// ��ÿ�����Ӷ���������
	void DivideArchiving();

	// ����ӵ�����б�
	void GetCrowd();
};

class GetGbest : MeshCrowd
{
public:
	double* archiveProbability;

	// Ĭ�Ϲ��캯��
	GetGbest() {}

	// ���캯��
	GetGbest(vector<Particle> curArchiveList, int meshDivCount, double* minList, double* maxList, int dim, int ParticleCount);

	// �������ӱ�ѡ��ĸ���
	void GetProbability();

	// 
	int GetGbestIndex();

	// ���㵱ǰ����Ⱥ��ÿһ�����ӵ�ȫ�����Ž⼯
	Particle* getGbest();



	// �������ӱ�ѡ��ĸ���
	void GetProbability1();

	// ������������ӣ�ӵ���ȸߵ����ӱ�����ĸ���Խ��
	vector<int> GetClearIndex(int threshN);

	// Clear����
	vector<Particle> Clear(int thresh);
};

