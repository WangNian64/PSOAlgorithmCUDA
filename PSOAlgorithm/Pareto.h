#pragma once
#include "PSO.h"
class Pareto
{
public:
	vector<Particle> particles;	//�������Ӽ���,���Ա���

	int cursor;				//��ǰ���ӵ��±꣬��ʼΪ-1
	int badSolutionNum;		//���Ž����
	int len;				//��������
public:
	// Ĭ�Ϲ��캯��
	Pareto() {}

	// ���캯��
	Pareto(vector<Particle> particles);

	// ��������
	~Pareto();

	//��ȡ��ǰ�±������
	Particle Next();

	//�ж��Ƿ��鿴����������
	bool HasNext();

	//ɾ���ӽ�
	void RemoveBadParticle();

	//�����������ӣ�����pareto���Ž⼯
	vector<Particle> GetPareto();

};
