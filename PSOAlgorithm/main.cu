#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>

#include "FitnessFunction.h"
int main()
{
	int problemSum = 1;//�������Ŀ
	int testSum = 1;//ÿ��ʵ���ܵ�ʵ�����
	for (int curProblem = 0; curProblem < problemSum; curProblem++)//�ܶ������
	{
		cout << "�ܵ�" + to_string(curProblem + 1) + "������" << endl;

#pragma region ����PSO����
		ifstream inputFile;
		inputFile.open("../../InputParas/InputPara" + to_string(curProblem + 1) + ".txt");

		//CPU
		ProblemParas proParas(inputFile);					// ��ʼ�������豸��ز��� 
		cout << "�ļ���ȡ�ɹ�" << endl;
		int ThreadsPerBlock = 100;							// һ��Block��100��thread
		int BlockSum = 1;									// Block����Ŀ
		int dim = proParas.DeviceSum * 3;					// ��ά��=�豸��*3(x,y,����)
		//CPU
		PSOPara psopara(dim);
		psopara.mesh_div_count = 4;							// ���񻮷���Ŀ
		psopara.problemParas = proParas;					// ��������Ĳ���
		psopara.particle_num_ = ThreadsPerBlock * BlockSum;	// ���Ӹ���
		psopara.max_iter_num_ = 400;						// ����������
		psopara.fitness_count_ = 2;							// ��Ӧ����Ŀ
		psopara.archive_max_count = 50;						// archive����������Ŀ
		psopara.SetDt(1.0);									// ʱ�䲽��
		psopara.SetWstart(0.9);								// ��ʼȨ��
		psopara.SetWend(0.4);								// ����Ȩ��
		psopara.SetC1(1.49445);								// ���ٶ�����1
		psopara.SetC2(1.49445);								// ���ٶ�����2
		psopara.SetLowBound(0, 0, DeviceDirect::Default);	// position��������Χ����

		psopara.blockSum = BlockSum;
		psopara.threadsPerBlock = ThreadsPerBlock;
		//��Ҫ���豸����ȡ�����ֵ��ֻ��ȡ��3.��
		psopara.SetUpBound(proParas.workShopLength, proParas.workShopWidth, DeviceDirect::Rotate270 + 1);// position��������Χ����
#pragma endregion

#pragma region ����PSO�㷨����������
//GPU
		PSOOptimizer psooptimizer(&psopara, proParas);//PSO�㷨����
		string curProblemFolderName = "Problem" + to_string(curProblem + 1);
		for (int curTest = 0; curTest < testSum; curTest++) {//ÿ�������ܶ��
			clock_t startTime, endTime;//��¼����ʱ��

			startTime = clock();//��ʱ��ʼ
#pragma region ��ʼ��
			psooptimizer.InitialAllParticles();//��ʼ���������� CPU
			psooptimizer.InitialArchiveList();//��ʼ��Archive�浵	CPU
			psooptimizer.InitGbest();//��ʼ��ȫ������		CPU
#pragma endregion

#pragma region ������������&��ÿһ�ε���Ӧ��ֵ
//Ŀ��1��ֵ����archiveList1�У�Ŀ��2��ֵ����archiveList2��
//��n��ʵ��ŵ��ļ���ȥ
//�ļ��е����ֽ�Testn
			ofstream OutFile;
			ofstream OutFile1;
			string curTestFolderName = "Test" + to_string(curTest + 1);
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList1.txt");
			OutFile1.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/archiveList2.txt");
			for (int i = 0; i < psooptimizer.max_iter_num_; i++)//��ʼ���в���
			{
				cout << (i + 1) << endl;
				psooptimizer.UpdateAllParticles();//�����������ӵ�λ�ú��ٶ�
				psooptimizer.UpdatePbest();//����pbest 

				psooptimizer.UpdateArchiveList();//�����ⲿ�浵����
				psooptimizer.UpdateGbest();//����gbest

				//�洢ÿ�ε�����Archive����
				double minFitness1, minFitness2;
				minFitness1 = minFitness2 = INT_MAX;
				cout << minFitness1 << endl;
				//archiveList��CPU��
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
			cout << "����" << psopara.max_iter_num_ << "�ε�������ʱ:" << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

#pragma region �����豸�ߴ�&���ղ��ֽ��&���ߵ������
			OutFile.open("../../Results/" + curProblemFolderName + "/" + curTestFolderName + "/FinalResult.txt");




#pragma region ��¼���ղ��ֽ��
			int resultIndex = 0;
			int minHandleCost = INT_MAX;
			int minConveyValue = INT_MAX;
			//����ѡ��������ɱ���͵�
			for (int i = 0; i < psooptimizer.archive_list.size(); i++)
			{
				if (psooptimizer.archive_list[i].fitness_[0] < minHandleCost)
				{
					minHandleCost = psooptimizer.archive_list[i].fitness_[0];
					resultIndex = i;
				}
			}
			//����ѡ���ͻ��ɱ��͵�
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

			//������һЩ������Ҫ��GPU->CPU
			//ֻ������Ҫ��
#pragma region ��¼�豸�ߴ�
			Vector2* deviceParaListSize_CPU = new Vector2[proParas.DeviceSum];
			//��ѭ���ķ�ʽ���������ݵ�CPU��ÿ��ֻ����һ��size�ĵ�ַ
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
					//�����������GPU�ϵ�
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
#pragma region ��¼��������꣨��ת֮��ģ������豸���꣩

			//����һ��inoutPoints
			//int ioPointsSize = psooptimizer.bestPathInfoList[fitnessIndex].inoutPSize;
			int ioPointsSize = psooptimizer.inoutPSize;//
			InoutPoint* ioPoints = new InoutPoint[ioPointsSize];
			cudaMemcpy(ioPoints, psooptimizer.curBestPath_InoutPoints, sizeof(InoutPoint) * ioPointsSize, cudaMemcpyDeviceToHost);

			OutFile << to_string(ioPointsSize) + "\n";//�������Ŀ
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

#pragma region ��¼�����·��
			////�ȴ�ÿ�ֻ����·������
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
			//	//����s2
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

#pragma region ��¼ֱ�����ͻ���ת�����ͻ�����
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

		//�����ڴ�

#pragma endregion
	}
	return 0;
}
