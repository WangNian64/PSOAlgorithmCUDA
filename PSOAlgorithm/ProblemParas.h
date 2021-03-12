#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "DevicePara.h"
using namespace std;
//�ַ����ָ��
static vector<string> split(const string& str, const string& pattern)
{
	vector<string> res;
	if ("" == str)
		return res;

	string strs = str + pattern;

	size_t pos = strs.find(pattern);
	size_t size = strs.size();

	while (pos != string::npos)
	{
		string x = strs.substr(0, pos);
		res.push_back(x);//stoi(x)ת����
		strs = strs.substr(pos + 1, size);
		pos = strs.find(pattern);
	}
	return res;
}


//�ɱ��������
struct CostPara {
	double MatFlow;		//��������
	double UnitCost;	//��λ�������ĳɱ�
};

struct ProblemParas
{
	//�ṹ��Ҫ�޸ģ�DevicePara��CargoType��Ҫ�ֽ�
	int DeviceSum;						//�豸����

	//DevicePara* deviceParaList;			//�豸�����б�
	int* ID;								//�豸ID
	double* workSpeed;						//�ӹ�/����1��λ���ϵ�ʱ��
	Vector2* size;							//�豸�ߴ磨�ֱ���x���y��ĳ��ȣ�
	Vector2* axis;							//�豸����
	DeviceDirect* direct;					//�豸����
	double* spaceLength;					//��϶��Ϊ��ʵ�־���Լ����

	//��adjPoint��Ŀ�йص�����
	int* adjPInCount;						//adjPointIn�������Ŀ����
	int* adjPOutCount;						//adjPointOut�������Ŀ����
	int* accumAdjPInCount;					//adjPointIn��Ŀ���ۼ�����
	int* accumAdjPOutCount;					//adjPointOut��Ŀ���ۼ�����

	//����ڵ�����飨��Ӱ�������ߵĲ��֣�
	int totalInPoint;						//��ڵ������Ŀ
	int totalOutPoint;						//���ڵ������Ŀ
	AdjPoint* adjPointsIn;	//���
	AdjPoint* adjPointsOut;	//����

	int horiPointCount;				//δȥ��ǰ����ˮƽ����ĵ����Ŀ
	int vertPointCount;				//δȥ��ǰ���д�ֱ����ĵ����Ŀ

	double workShopLength;				//���䳤��
	double workShopWidth;				//������

	Vector2 entrancePos;				//�ֿ��������	
	Vector2 exitPos;					//�ֿ��������

	//���ϲ����б�
	int CargoTypeNum;					//����������Ŀ
	//CargoType* cargoTypeList;			//���������б�
	int cargoTypeSize;					//���������б���ڴ��С
	int totalLinkSum;					//�ܵ���������Ŀ

	int* deviceSum;							//�������豸��Ŀ
	int* linkSum;							//�豸��Ե���Ŀ
	int* accumLinkSum;						//linkSum���ۼ���Ŀ
	DeviceLink* deviceLinkList;				//�豸�����б�
	double* totalVolume;					//�����ϵ���������

	int fixedLinkPointSum;					//ÿһ��link�Ĺ̶�����ĿΪ50(û��ȥ��֮ǰ��)
	int fixedUniqueLinkPointSum;			//ȥ�غ��ÿһ��link�Ĺ̶�����ĿΪ20
	//���ͻ�����
	double convey2DeviceDist;//���ͻ����豸�ľ��루Ѱ·��ʱ��Ҫ���ǣ�
	double conveyWidth;//���ͻ����
	double conveyMinLength;//���ͻ���̳���
	double conveySpeed;//���ͻ������ٶ�
	double strConveyorUnitCost;//��λֱ�����ͻ��ɱ�
	double curveConveyorUnitCost;//����ת�����ͻ��ɱ�

	double conveyMinDist;//������������֮�����̾���
	//��ǰ�Ĳ���
	//CostPara** costParaArray;		//�ɱ������������(�����������������ɱ�)
	//double** minDistArray;			//�豸��С��������
	//DeviceRelation** deviceRelateArray;//�豸���������

	ProblemParas() {}
	ProblemParas(ifstream& inputFile)
	{
		//��ȡ����
		string line;
		cout << "3333" << endl;
		if (inputFile) // �и��ļ�
		{
			horiPointCount = 0;
			vertPointCount = 0;

			fixedLinkPointSum = 50;
			fixedUniqueLinkPointSum = 20;

#pragma region �⼸������Ŀǰ����
			//costParaArray = new CostPara * [DeviceSum];
			//for (int i = 0; i < DeviceSum; i++) {
			//	costParaArray[i] = new CostPara[DeviceSum];
			//}

			//minDistArray = new double* [DeviceSum];
			//for (int i = 0; i < DeviceSum; i++) {
			//	minDistArray[i] = new double[DeviceSum];
			//}

			//deviceRelateArray = new DeviceRelation * [DeviceSum];
			//for (int i = 0; i < DeviceSum; i++) {
			//	deviceRelateArray[i] = new DeviceRelation[DeviceSum];
			//}
#pragma endregion

#pragma region ����ߴ�	
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			vector<string> shopSizeStr = split(line, ",");
			workShopLength = atof(shopSizeStr[0].c_str());
			workShopWidth = atof(shopSizeStr[1].c_str());
#pragma endregion

#pragma region �ֿ��������	
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			vector<string> enterPosStr = split(line, ",");
			entrancePos.x = atof(enterPosStr[0].c_str());
			entrancePos.y = atof(enterPosStr[1].c_str());
#pragma endregion

#pragma region �ֿ��������	
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			vector<string> exitPosStr = split(line, ",");
			exitPos.x = atof(exitPosStr[0].c_str());
			exitPos.y = atof(exitPosStr[1].c_str());
#pragma endregion

#pragma region �豸��ز���
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			DeviceSum = atoi(line.c_str());//�豸��Ŀ

			//deviceParaList = new DevicePara[DeviceSum];
			ID = new int[DeviceSum]; 				//�豸ID
			workSpeed = new double[DeviceSum];		//�ӹ�/����1��λ���ϵ�ʱ��
			size = new Vector2[DeviceSum];			//�豸�ߴ磨�ֱ���x���y��ĳ��ȣ�
			axis = new Vector2[DeviceSum];			//�豸����
			direct = new DeviceDirect[DeviceSum];	//�豸����
			spaceLength = new double[DeviceSum];	//��϶��Ϊ��ʵ�־���Լ����

			//��adjPoint��Ŀ�йص�����
			adjPInCount = new int[DeviceSum];						//adjPointIn�������Ŀ����
			adjPOutCount = new int[DeviceSum];						//adjPointOut�������Ŀ����

			accumAdjPInCount = new int[DeviceSum];					//adjPointIn��Ŀ���ۼ�����
			accumAdjPOutCount = new int[DeviceSum];					//adjPointOut��Ŀ���ۼ�����

			totalInPoint = totalOutPoint = 0;//��Ҫ�ۼ�

			//��Ҫ�ȱ���һ�飬Ȼ�������¸�ֵһ�β���
			//����vector����Щ����
			vector<AdjPoint> tempAdjPInList;
			vector<AdjPoint> tempAdjPOutList;
			for (int i = 0; i < DeviceSum; i++)
			{
				//int curAdjInCount = 0; //��һ�ֵ�AdjIn����Ŀ
				//int curAdjOutCount = 0;//��һ�ֵ�AdjOut����Ŀ 
				//��temp vector��ת
				//vector<AdjPoint> tempAdjPInList;
				//vector<AdjPoint> tempAdjPOutList;

				getline(inputFile, line);
				vector<string> strSplit = split(line, " ");

				ID[i] = atoi(strSplit[0].c_str());
				workSpeed[i] = atof(strSplit[1].c_str());
				vector<string> sizeStr = split(strSplit[2], ",");
				size[i].x = atof(sizeStr[0].c_str());
				size[i].y = atof(sizeStr[1].c_str());

				spaceLength[i] = atof(strSplit[3].c_str());
				//�ڽӵ��λ��
				if (strSplit[4] != "null")
				{
					vector<string> adjStrList = split(strSplit[4], "|");//�����ÿ���㣬�����ÿ��in��out�ܵĵ㣬����Ҫ�ֿ�
					int InIndex, OutIndex;
					InIndex = OutIndex = 0;
					for (int k = 0; k < adjStrList.size(); k++)
					{
						vector<string> adjStr = split(adjStrList[k], ",");//�ٴη���
						AdjPoint adjPoint;
						adjPoint.direct = adjPoint.GetDirect(adjStr[0]);
						adjPoint.inoutType = (InoutType)(atoi(adjStr[1].c_str()) - 1);
						adjPoint.pos.x = atof(adjStr[2].c_str());
						adjPoint.pos.y = atof(adjStr[3].c_str());
						if (adjPoint.inoutType == 0) {//����ڵ�
							adjPoint.index = InIndex++;
							tempAdjPInList.push_back(adjPoint);
						}
						else {//�ǳ��ڵ�
							adjPoint.index = OutIndex++;
							tempAdjPOutList.push_back(adjPoint);
						}
					}
					adjPInCount[i] = InIndex;
					adjPOutCount[i] = OutIndex;

					accumAdjPInCount[i] = totalInPoint;//�ۼӵĵ�����ô�����(��һ����0)
					accumAdjPOutCount[i] = totalOutPoint;
					//�ۼ�in��outPoint��Ŀ
					totalInPoint = totalInPoint + adjPInCount[i];
					totalOutPoint = totalOutPoint + adjPOutCount[i];
				}
			}
			//����ڵ�����飨��Ӱ�������ߵĲ��֣�
			adjPointsIn = new AdjPoint[totalInPoint];	//���
			cout << totalInPoint << endl;
			adjPointsOut = new AdjPoint[totalOutPoint];	//����
			cout << totalOutPoint << endl;

			//��vector�ĵ㸴�Ƶ�*����
			for (int i = 0; i < tempAdjPInList.size(); i++)
			{
				adjPointsIn[i] = tempAdjPInList[i];
			}
			for (int i = 0; i < tempAdjPOutList.size(); i++)
			{
				adjPointsOut[i] = tempAdjPOutList[i];
			}

#pragma region ����ˮƽ�ʹ�ֱ�����Ŀ��δ�����غϵĵ㣩
			//���ˮƽ�ʹ�ֱ����ڵ����Ŀ��һ����horiCount��vertCount��
			for (int i = 0; i < DeviceSum; i++)
			{
				for (int pointIndex = 0; pointIndex < adjPInCount[i]; pointIndex++)
				{
					AdjPoint& p = adjPointsIn[accumAdjPInCount[i] + pointIndex];//���±�=ǰ���ۼӵ�+index
					if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//����
					{
						horiPointCount++;
					}
					else {//����
						vertPointCount++;
					}
				}
				for (int pointIndex = 0; pointIndex < adjPOutCount[i]; pointIndex++)
				{
					AdjPoint& p = adjPointsOut[accumAdjPOutCount[i] + pointIndex];//���±�=ǰ���ۼӵ�+index
					if (p.direct == PointDirect::Up || p.direct == PointDirect::Down)//����
					{
						horiPointCount++;
					}
					else {//����
						vertPointCount++;
					}
				}
			}
			//�ֿ����&����2���㣨horiCount��vertCount+=2��
			horiPointCount += 2;
			vertPointCount += 2;
			//ÿ���豸������ĸ���Χ
			horiPointCount = horiPointCount + 4 * DeviceSum;
			vertPointCount = vertPointCount + 4 * DeviceSum;
#pragma endregion

#pragma endregion

#pragma region �����߲���
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			vector<string> conveyInfoStr = split(line, " ");
			conveyWidth = atof(conveyInfoStr[0].c_str());
			conveySpeed = atof(conveyInfoStr[1].c_str());
			conveyMinLength = atof(conveyInfoStr[2].c_str());
			convey2DeviceDist = atof(conveyInfoStr[3].c_str());
			strConveyorUnitCost = atof(conveyInfoStr[4].c_str());
			curveConveyorUnitCost = atof(conveyInfoStr[5].c_str());
			conveyMinDist = conveyMinLength + conveyWidth;
#pragma endregion

#pragma region ���ϲ���
			getline(inputFile, line);//��һ��
			getline(inputFile, line);
			CargoTypeNum = atoi(line.c_str());//����������Ŀ

			cout << CargoTypeNum << endl;

			//cargoTypeList = new CargoType[CargoTypeNum];
			deviceSum = new int[CargoTypeNum];							//�������豸��Ŀ
			linkSum = new int[CargoTypeNum];							//�豸��Ե���Ŀ
			accumLinkSum = new int[CargoTypeNum];						//linkSum���ۼ���Ŀ
			totalVolume = new double[CargoTypeNum];					//�����ϵ���������

			vector<CargoType> tempCargoTypeList;
			totalLinkSum = 0;


			cout << "111" << endl;
			vector<DeviceLink> tempDeviceLinkList;
			for (int i = 0; i < CargoTypeNum; i++)
			{

				getline(inputFile, line);
				vector<string> strSplit = split(line, " ");
				linkSum[i] = atoi(strSplit[0].c_str());
				deviceSum[i] = linkSum[i] + 1;

				accumLinkSum[i] = totalLinkSum;
				//deviceLinkList = new DeviceLink[cargoTypeList[i].linkSum];

				for (int j = 0; j < linkSum[i]; j++)
				{
					DeviceLink tempLink = DeviceLink();
					vector<string> deviceLinkStr = split(strSplit[j + 1], "-");
					//��Ϊ��ںͳ���
					if (deviceLinkStr[0] == "ENTER")//˵���ǲֿ����
					{
						tempLink.outDeviceIndex = -1;
					}
					else
					{
						vector<string> devicePointStr = split(deviceLinkStr[0], ",");
						tempLink.outDeviceIndex = atoi(devicePointStr[0].c_str()) - 1;
						tempLink.outPointIndex = atoi(devicePointStr[1].c_str()) - 1;
					}
					if (deviceLinkStr[1] == "EXIT")//˵���ǲֿ����
					{
						tempLink.inDeviceIndex = -2;
					}
					else
					{
						vector<string> devicePointStr = split(deviceLinkStr[1], ",");
						tempLink.inDeviceIndex = atoi(devicePointStr[0].c_str()) - 1;
						tempLink.inPointIndex = atoi(devicePointStr[1].c_str()) - 1;
					}
					tempDeviceLinkList.push_back(tempLink);
				}
				totalVolume[i] = atof(strSplit[strSplit.size() - 1].c_str());
			}
			totalLinkSum = tempDeviceLinkList.size();
			cout << totalLinkSum << endl;
			deviceLinkList = new DeviceLink[totalLinkSum];
			for (int i = 0; i < tempDeviceLinkList.size(); i++)
			{
				deviceLinkList[i] = tempDeviceLinkList[i];
			}
			//vector������*����
#pragma endregion


			cout << "111" << endl;


#pragma region ��λ����������ɱ�����
			//getline(inputFile, line);//��һ��
			//for (int i = 0; i < DeviceSum; i++)
			//{
			//	getline(inputFile, line);
			//	vector<string> strSplit = split(line, ",");
			//	for (int j = 0; j < DeviceSum; j++)
			//	{
			//		costParaArray[i][j].UnitCost = atof(strSplit[j].c_str());
			//	}
			//}
#pragma endregion

#pragma region ����������
//getline(inputFile, line);//��һ��
//for (int i = 0; i < DeviceSum; i++) {
//	getline(inputFile, line);
//	vector<string> strSplit = split(line, ",");
//	for (int j = 0; j < DeviceSum; j++)
//	{
//		costParaArray[i][j].MatFlow = atof(strSplit[j].c_str());
//	}
//}
#pragma endregion

#pragma region �豸��С�������� 
//getline(inputFile, line);//��һ��
//for (int i = 0; i < DeviceSum; i++) {
//	getline(inputFile, line);
//	vector<string> strSplit = split(line, ",");
//	for (int j = 0; j < DeviceSum; j++)
//	{
//		minDistArray[i][j] = atof(strSplit[j].c_str());
//	}
//}
#pragma endregion

#pragma region �豸�������������
//getline(inputFile, line);//��һ��
//for (int i = 0; i < DeviceSum; i++) {
//	getline(inputFile, line);
//	vector<string> strSplit = split(line, ",");
//	for (int j = 0; j < DeviceSum; j++)
//	{
//		deviceRelateArray[i][j] = (DeviceRelation)atoi(strSplit[j].c_str());
//	}
//}
#pragma endregion

		}
		else // û�и��ļ�
		{
			cout << "no such file" << endl;
		}

#pragma region �����豸ͼ�ṹ(�ƺ�û��ʲô��)
		////�㼯
		//for (int i = 0; i < DeviceSum; i++)
		//{
		//	deviceGraph.InsertVertex(i);
		//}
		////��
		//for (int i = 0; i < CargoTypeNum; i++)
		//{
		//	for (int j = 0; j < cargoTypeList[i].deviceSum - 1; j++)
		//	{
		//		int firstEdge = cargoTypeList[i].deviceList[j] - 1;//�豸���-1
		//		int secondEdge = cargoTypeList[i].deviceList[j + 1] - 1;
		//		deviceGraph.InsertEdge(firstEdge, secondEdge);
		//	}
		//}
		//deviceGraph.ShowGraph();
#pragma endregion
		cout << "��ȡ���" << endl;

		cout << "111" << endl;
	}

	//ProblemParas(const ProblemParas & para)
	//{
	//	this->DeviceSum = para.DeviceSum;
	//	this->deviceParaList = new DevicePara[this->DeviceSum];
	//	for (int i = 0; i < DeviceSum; i++)
	//	{
	//		this->deviceParaList[i] = para.deviceParaList[i];
	//	}
	//	//����ߴ�
	//	this->workShopLength = para.workShopLength;
	//  this->workShopWidth = para.workShopWidth;
	//	//�ֿ��������
	//	this->entrancePos = para.entrancePos;
	//	this->costParaArray = new CostPara * [this->DeviceSum];
	//	this->minDistArray = new double * [this->DeviceSum];
	//	this->deviceRelateArray = new DeviceRelation * [this->DeviceSum];
	//	for (int i = 0; i < this->DeviceSum; i++)
	//	{
	//		this->costParaArray[i] = new CostPara[this->DeviceSum];
	//		this->minDistArray[i] = new double[this->DeviceSum];
	//		this->deviceRelateArray[i] = new DeviceRelation[this->DeviceSum];
	//	}
	//	for (int i = 0; i < this->DeviceSum; i++)
	//	{
	//		for (int j = 0; j < this->DeviceSum; j++)
	//		{
	//			this->costParaArray[i][j] = para.costParaArray[i][j];
	//			this->minDistArray[i][j] = para.minDistArray[i][j];
	//			this->deviceRelateArray[i][j] = para.deviceRelateArray[i][j];
	//		}
	//	}
	//	//���ϲ����б�
	//	this->CargoTypeNum = para.CargoTypeNum;
	//	this->cargoTypeList = new CargoType[this->CargoTypeNum];

	//	this->deviceGraph = para.deviceGraph;//�豸����ͼ�ṹ

	//	this->conveySpeed = para.conveySpeed;//���ͻ������ٶ�
	//	this->spaceLength = para.spaceLength;//���豸�������һ��ǣ�Ϊ�˸��õľ���Լ����
	//}
};
