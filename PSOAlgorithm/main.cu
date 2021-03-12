#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>

#include "FitnessFunction.h"
int main()
{
	int problemSum = 1;//问题的数目
	int testSum = 1;//每个实验跑的实验次数
	for (int curProblem = 0; curProblem < problemSum; curProblem++)//跑多个问题
	{
		cout << "跑第" + to_string(curProblem + 1) + "个问题" << endl;

#pragma region 设置PSO参数
		ifstream inputFile;
		inputFile.open("../../InputParas/InputPara" + to_string(curProblem + 1) + ".txt");

		//CPU
		ProblemParas proParas(inputFile);					// 初始化所有设备相关参数 
		cout << "文件读取成功" << endl;
		int ThreadsPerBlock = 100;							// 一个Block中100个thread
		int BlockSum = 1;									// Block的数目
		int dim = proParas.DeviceSum * 3;					// 总维度=设备数*3(x,y,朝向)
		//CPU
		PSOPara psopara(dim);
		psopara.mesh_div_count = 4;							// 网格划分数目
		psopara.problemParas = proParas;					// 布局问题的参数
		psopara.particle_num_ = ThreadsPerBlock * BlockSum;	// 粒子个数
		psopara.max_iter_num_ = 400;						// 最大迭代次数
		psopara.fitness_count_ = 2;							// 适应度数目
		psopara.archive_max_count = 50;						// archive数组的最大数目
		psopara.SetDt(1.0);									// 时间步长
		psopara.SetWstart(0.9);								// 初始权重
		psopara.SetWend(0.4);								// 结束权重
		psopara.SetC1(1.49445);								// 加速度因子1
		psopara.SetC2(1.49445);								// 加速度因子2
		psopara.SetLowBound(0, 0, DeviceDirect::Default);	// position的搜索范围下限

		psopara.blockSum = BlockSum;
		psopara.threadsPerBlock = ThreadsPerBlock;
		//不要让设备朝向取到最大值，只能取到3.几
		psopara.SetUpBound(proParas.workShopLength, proParas.workShopWidth, DeviceDirect::Rotate270 + 1);// position的搜索范围上限
#pragma endregion

#pragma region 调用PSO算法，并输出结果
//GPU
		PSOOptimizer psooptimizer(&psopara, proParas);//PSO算法对象
		string curProblemFolderName = "Problem" + to_string(curProblem + 1);
		for (int curTest = 0; curTest < testSum; curTest++) {//每个问题跑多次
			clock_t startTime, endTime;//记录调用时间

			startTime = clock();//计时开始
#pragma region 初始化
			psooptimizer.InitialAllParticles();//初始化所有粒子 CPU
			psooptimizer.InitialArchiveList();//初始化Archive存档	CPU
			psooptimizer.InitGbest();//初始化全局最优		CPU
#pragma endregion

#pragma region 迭代更新粒子&存每一次的适应度值
//目标1的值放在archiveList1中，目标2的值放在archiveList2中
//第n次实验放到文件里去
//文件夹的名字叫Testn
			ofstream OutFile;
			ofstream OutFile1;
			string curTestFolderName = "Test" + to_string(curTest + 1);
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList1.txt");
			OutFile1.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList2.txt");
			for (int i = 0; i < psooptimizer.max_iter_num_; i++)//开始并行操作
			{
				cout << (i + 1) << endl;
				psooptimizer.UpdateAllParticles();//更新所有粒子的位置和速度
				psooptimizer.UpdatePbest();//更新pbest 

				psooptimizer.UpdateArchiveList();//更新外部存档集合
				psooptimizer.UpdateGbest();//更新gbest

				//存储每次迭代的Archive集合
				double minFitness1, minFitness2;
				minFitness1 = minFitness2 = INT_MAX;
				cout << minFitness1 << endl;
				//archiveList在CPU端
				for (auto it = psooptimizer.archive_list.begin(); it != psooptimizer.archive_list.end(); it++)
				{
					minFitness1 = min(minFitness1, it->fitness_[0]);
				}
				string f1line = to_string(minFitness1) + "\n";
				OutFile << f1line;
				cout << f1line << endl;

				for (auto it = psooptimizer.archive_list.begin(); it != psooptimizer.archive_list.end(); it++)
				{
					minFitness2 = min(minFitness2, it->fitness_[1]);
				}
				string f2line = to_string(minFitness2) + "\n";
				OutFile1 << f2line;
				cout << f2line << endl;
			}
			OutFile.close();
			OutFile1.close();
#pragma endregion

			endTime = clock();
			cout << "迭代" << psopara.max_iter_num_ << "次的最终用时:" << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

#pragma region 保存设备尺寸&最终布局结果&连线点的坐标
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/FinalResult.txt");




#pragma region 记录最终布局结果
			int resultIndex = 0;
			int minHandleCost = INT_MAX;
			int minConveyValue = INT_MAX;
			//优先选物料运输成本最低的
			for (int i = 0; i < psooptimizer.archive_list.size(); i++)
			{
				if (psooptimizer.archive_list[i].fitness_[0] < minHandleCost)
				{
					minHandleCost = psooptimizer.archive_list[i].fitness_[0];
					resultIndex = i;
				}
			}
			//优先选输送机成本低的
			//for (int i = 0; i < psooptimizer.archive_list.size(); i++)
			//{
			//	if (psooptimizer.archive_list[i].fitness_[1] < minConveyValue)
			//	{
			//		minConveyValue = psooptimizer.archive_list[i].fitness_[0];
			//		resultIndex = i;
			//	}
			//}

			for (int i = 0; i < dim; i += 3)
			{
				OutFile << psooptimizer.archive_list[resultIndex].position_[i];
				OutFile << ",";
				OutFile << psooptimizer.archive_list[resultIndex].position_[i + 1];
				OutFile << "\n";
			}
#pragma endregion

			//后面有一些参数需要从GPU->CPU
			//只复制需要的
#pragma region 记录设备尺寸
			Vector2* deviceParaListSize_CPU = new Vector2[proParas.DeviceSum];
			//用循环的方式来复制数据到CPU，每次只复制一个size的地址
			for (int i = 0; i < proParas.DeviceSum; i++)
			{
				cudaMemcpy(deviceParaListSize_CPU + i, &psooptimizer.problemParas.size[i], sizeof(Vector2), cudaMemcpyDeviceToHost);
			}

			for (int i = 2; i < dim; i += 3)
			{
				DeviceDirect direct = (DeviceDirect)(int)psooptimizer.archive_list[resultIndex].position_[i];
				string line = "";
				if (direct == DeviceDirect::Rotate90 || direct == DeviceDirect::Rotate270)
				{
					//这里的数据是GPU上的
					line = to_string(deviceParaListSize_CPU[i / 3].y) + "," +
						to_string(deviceParaListSize_CPU[i / 3].x);
				}
				else {
					line = to_string(deviceParaListSize_CPU[i / 3].x) + "," +
						to_string(deviceParaListSize_CPU[i / 3].y);
				}
				OutFile << line + "\n";
			}
#pragma endregion

			int fitnessIndex = 0;
#pragma region 记录出入口坐标（旋转之后的，不带设备坐标）

			//复制一遍inoutPoints
			//int ioPointsSize = psooptimizer.bestPathInfoList[fitnessIndex].inoutPSize;
			int ioPointsSize = psooptimizer.inoutPSize;//
			InoutPoint* ioPoints = new InoutPoint[ioPointsSize];
			cudaMemcpy(ioPoints, psooptimizer.curBestPath_InoutPoints, sizeof(InoutPoint) * ioPointsSize, cudaMemcpyDeviceToHost);

			OutFile << to_string(ioPointsSize) + "\n";//出入口数目
			for (int i = 0; i < ioPointsSize; i++)
			{
				if (ioPoints[i].pointDirect == PointDirect::Up || ioPoints[i].pointDirect == PointDirect::Down)
				{
					OutFile << "Vertical ";
				}
				else
				{
					OutFile << "Horizon ";
				}
				OutFile << ioPoints[i].pointAxis.x;
				OutFile << " ";
				OutFile << ioPoints[i].pointAxis.y;
				OutFile << "\n";
			}
#pragma endregion

#pragma region 记录出入口路径
			////先存每种货物的路径条数
			//string line = "";
			//for (int i = 0; i < proParas.CargoTypeNum; i++)
			//{
			//	line += to_string(proParas.cargoTypeList[i].deviceSum - 1);
			//	if (i != proParas.CargoTypeNum - 1)
			//	{
			//		line += " ";
			//	}
			//}
			//OutFile << line << "\n";

			//vector<PointLink> p = psooptimizer.archive_list[resultIndex].pointLinks;
			//for (int i = 0; i < p.size(); i++)
			//{
			//	string s1, s2;
			//	DevicePara device1, device2;

			//	//s1 = to_string(p[i].device1Index) + " " + to_string(p[i].device2Index);
			//	OutFile << to_string(p[i].device1Index) + " " + to_string(p[i].device2Index) + " ";
			//	//计算s2
			//	for (int j = 0; j < p[i].points.size(); j++)
			//	{
			//		OutFile /*<< fixed << setprecision(1)*/ << p[i].points[j].x;
			//		OutFile << ",";
			//		OutFile /*<< fixed << setprecision(1)*/ << p[i].points[j].y;
			//		//s2 += to_string(p[i].points[j].x) + "," + to_string(p[i].points[j].y);
			//		if (j != p[i].points.size() - 1)
			//		{
			//			//s2 += "|";
			//			OutFile << "|";
			//		}
			//	}
			//	OutFile << "\n";
			//	//string line = s1 + " " + s2 + "\n";
			//	//OutFile << line;
			//}
#pragma endregion

#pragma region 记录直线输送机和转弯输送机参数
//GPU->CPU
			int strInfoListSum = psooptimizer.curBestPath_StrConveyorListSum[0];
			StraightConveyorInfo* strInfoList = new StraightConveyorInfo[strInfoListSum];
			cudaMemcpy(strInfoList, psooptimizer.curBestPath_StrConveyorList, sizeof(StraightConveyorInfo) * strInfoListSum, cudaMemcpyDeviceToHost);

			int curveInfoListSum = psooptimizer.curBestPath_CurveConveyorListSum[0];
			Vector2Int* curveInfoList = new Vector2Int[curveInfoListSum];
			cudaMemcpy(curveInfoList, psooptimizer.curBestPath_CurveConveyorList, sizeof(Vector2Int)* curveInfoListSum, cudaMemcpyDeviceToHost);

			OutFile << strInfoListSum << "\n";
			for (int i = 0; i < strInfoListSum; i++)
			{
				OutFile << to_string(strInfoList[i].startPos.x) << "," << to_string(strInfoList[i].startPos.y)
					<< ";" << to_string(strInfoList[i].endPos.x) << "," << to_string(strInfoList[i].endPos.y)
					<< ";" << to_string(strInfoList[i].startHnum) << ";" << to_string(strInfoList[i].startVnum)
					<< ";" << to_string(strInfoList[i].endHnum) << ";" << to_string(strInfoList[i].endVnum)
					<< "\n";
			}
			OutFile << curveInfoListSum << "\n";
			for (int i = 0; i < curveInfoListSum; i++)
			{
				OutFile << to_string(curveInfoList[i].x) << "," << to_string(curveInfoList[i].y) << "\n";
			}
#pragma endregion

			OutFile.close();
#pragma endregion

		}

		//回收内存

#pragma endregion
	}
	return 0;
}
