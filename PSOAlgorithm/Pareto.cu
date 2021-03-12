#pragma once
#include "Pareto.h"
static bool Judge_ParetoOptimal(Particle curParticle, vector<Particle> particles, int cursor);

static bool CompareFitness(Particle curParticle, Particle otherParticle);
// ���캯��
Pareto::Pareto(vector<Particle> particles)
{
	this->particles = particles;
	this->cursor = -1;
	this->len = particles.size();
	this->badSolutionNum = 0;
}

// ��������
Pareto::~Pareto()
{

}

//��ȡ��ǰ�±������
Particle Pareto::Next()
{
	this->cursor = this->cursor + 1;
	return particles[this->cursor];
}

//�ж��Ƿ��鿴����������
bool Pareto::HasNext()
{
	return len > cursor + 1 + badSolutionNum;
}

//ɾ����ǰ�±���ӽ�
void Pareto::RemoveBadParticle()
{
	//ɾ���ӽ�
	particles.erase(particles.begin() + cursor);
	cursor--;
	badSolutionNum++;
}

//�����������ӣ�����pareto���Ž⼯
vector<Particle> Pareto::GetPareto()
{
	while (HasNext())
	{
		//��ȡ��ǰ����
		Particle curParticle = Next();
		//�жϵ�ǰ�����Ƿ�pareto����
		if (Judge_ParetoOptimal(curParticle, particles, cursor) == false)
		{
			//cout << cursor << endl;
			RemoveBadParticle();
		}
	}
	return particles;
}

//�ж������Ƿ�Ϊpareto���Ž⣨�������������ӶԱȣ�
static bool Judge_ParetoOptimal(Particle curParticle, vector<Particle> particles, int cursor)
{
	for (int i = 0; i < particles.size(); i++)
	{
		if (i == cursor)
			continue;
		//������ݼ��д���һ���������ȫ֧�䵱ǰ�⣬��ǰ��Ϊ�ӽ⣨�Ա���Ӧ�ȣ�
		if (CompareFitness(curParticle, particles[i]) == false)
		{
			return false;
		}
	}
	return true;
}

//�жϵ�ǰ���Ƿ��ܱ���ȫ֧��,false����ȫ֧��
static bool CompareFitness(Particle curParticle, Particle otherParticle)
{
	for (int i = 0; i < otherParticle.fitnessCount; i++)
	{
		//��Ӧ��ԽСԽ�ã���ȻҲ�����޸�
		if (curParticle.fitness_[i] < otherParticle.fitness_[i])
		{
			return true;
		}
	}
	return false;
}