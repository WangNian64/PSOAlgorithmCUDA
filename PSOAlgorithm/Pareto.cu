#pragma once
#include "Pareto.h"
static bool Judge_ParetoOptimal(Particle curParticle, vector<Particle> particles, int cursor);

static bool CompareFitness(Particle curParticle, Particle otherParticle);
// 构造函数
Pareto::Pareto(vector<Particle> particles)
{
	this->particles = particles;
	this->cursor = -1;
	this->len = particles.size();
	this->badSolutionNum = 0;
}

// 析构函数
Pareto::~Pareto()
{

}

//获取当前下标的粒子
Particle Pareto::Next()
{
	this->cursor = this->cursor + 1;
	return particles[this->cursor];
}

//判断是否检查看了所有粒子
bool Pareto::HasNext()
{
	return len > cursor + 1 + badSolutionNum;
}

//删除当前下标的劣解
void Pareto::RemoveBadParticle()
{
	//删除劣解
	particles.erase(particles.begin() + cursor);
	cursor--;
	badSolutionNum++;
}

//遍历所有粒子，更新pareto最优解集
vector<Particle> Pareto::GetPareto()
{
	while (HasNext())
	{
		//获取当前粒子
		Particle curParticle = Next();
		//判断当前粒子是否pareto最优
		if (Judge_ParetoOptimal(curParticle, particles, cursor) == false)
		{
			//cout << cursor << endl;
			RemoveBadParticle();
		}
	}
	return particles;
}

//判断粒子是否为pareto最优解（和其他所有粒子对比）
static bool Judge_ParetoOptimal(Particle curParticle, vector<Particle> particles, int cursor)
{
	for (int i = 0; i < particles.size(); i++)
	{
		if (i == cursor)
			continue;
		//如果数据集中存在一个解可以完全支配当前解，则当前解为劣解（对比适应度）
		if (CompareFitness(curParticle, particles[i]) == false)
		{
			return false;
		}
	}
	return true;
}

//判断当前解是否能被完全支配,false是完全支配
static bool CompareFitness(Particle curParticle, Particle otherParticle)
{
	for (int i = 0; i < otherParticle.fitnessCount; i++)
	{
		//适应度越小越好，当然也可以修改
		if (curParticle.fitness_[i] < otherParticle.fitness_[i])
		{
			return true;
		}
	}
	return false;
}