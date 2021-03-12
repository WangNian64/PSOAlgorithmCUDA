#pragma once
#include "PSO.h"
class Pareto
{
public:
	vector<Particle> particles;	//所有粒子集合,可以保留

	int cursor;				//当前粒子的下标，初始为-1
	int badSolutionNum;		//非优解个数
	int len;				//粒子总数
public:
	// 默认构造函数
	Pareto() {}

	// 构造函数
	Pareto(vector<Particle> particles);

	// 析构函数
	~Pareto();

	//获取当前下标的粒子
	Particle Next();

	//判断是否检查看了所有粒子
	bool HasNext();

	//删除劣解
	void RemoveBadParticle();

	//遍历所有粒子，更新pareto最优解集
	vector<Particle> GetPareto();

};
