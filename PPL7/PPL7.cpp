#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cfloat>
#include <streambuf>
#include <string>
#include <iomanip>
#include <chrono>
#include <vector>
#include <thread>

using namespace std;

/*
typedef struct
{
	float x;
	float y;
} Point;*/

typedef struct
{
	int v;
	float m;
	int pj;
} minResp;

const int OPER_UPDATE = 1;
const int OPER_BUILD = 2;
const int OPER_MIN = 3;
const int OPER_EXIT = 4;
const int OPER_INIT = 5;

void readGraph();
void writeMST();
void initWorkers();
void release();
void buildMST();
void runBuildKey(bool* key, bool* mstSet);
void runMinKey(bool* key, bool* mstSet);
void buildKey(bool* key, bool* mstSet, int pfrom, int count);
void minKey(int& u, float& m, int& pj, bool* key, bool* mstSet, int pfrom, int count);
//void writeLog();

int procCount = 3;
int curRank;
int n; // vertex count
int pathCount;
float* graph;
int* result;
int* minKeys;
int* jj;
float* mm;
//Point* points;
bool* key;
bool* mstSet;
int* partsSize;
int* starts;
int partFrom;
int partSize;

void worker();

int main(int argc, char* argv[])
{
	int ierr;
	ierr = MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) return ierr;
	MPI_Comm_rank(MPI_COMM_WORLD, &curRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	if (curRank == 0) {
		readGraph();
		minKeys = new int[procCount];
		mm = new float[procCount];
		jj = new int[procCount];
		key = new bool[n];
		mstSet = new bool[n];
		result = new int[n];
		initWorkers();
		buildMST();
		writeMST();
		release();
	}
	else {
		worker();
	}
	MPI_Finalize();
}

void initWorkers() {
	int nominalCount = n / (procCount - 1);
	partsSize = new int[procCount];
	starts = new int[procCount];
	for (int startIndex = 0, i = 1; startIndex < n; startIndex += nominalCount, i++) {
		int count = nominalCount;
		if (startIndex + nominalCount < n && startIndex + nominalCount * 2 >= n) {
			count = n - startIndex;
		}
		partsSize[i] = count;
		starts[i]] = startIndex;
	}
	for (int i = 1; i < procCount; i++) {
		int oper = OPER_INIT;
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&starts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&partsSize[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(graph, n * n, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
		key = new bool[n];
		mstSet = new bool[n];
	}
}

void worker() {
	int oper;
	MPI_Status status;
	MPI_Recv(&oper, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	switch (oper) {
	case OPER_INIT:
		MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&partFrom, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&partSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		graph = new float[n * n];
		MPI_Recv(graph, n * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		break;
	case OPER_EXIT:
		return;
		break;
	case OPER_BUILD:
		break;
	case OPER_MIN:
		break;
	case OPER_UPDATE:
		break;
	}
}

void readGraph()
{
	cout << "Reading graph ..." << endl;
	ifstream inf("graph.bin", ios::in | ios::binary);
	if (!inf)
	{
		cerr << "Graph reading error!" << endl;
		exit(1);
	}
	inf.read((char*)&n, sizeof(n));
	pathCount = (n * (n - 1)) / 2;
	float* paths = new float[pathCount];
	inf.read((char*)paths, pathCount * sizeof(float));
	//points = new Point[n];
	//inf.read((char*)points, n * sizeof(Point));
	inf.close();
	graph = new float[n * n];
	int p = 0;
	for (int i = 0; i < n; i++)
	{
		int j = 0;
		for (; j < i; j++)
		{
			graph[i * n + j] = graph[j * n + i];
		}
		j = i;
		graph[i * n + j] = 0;
		for (j++; j < n; j++)
		{
			graph[i * n + j] = paths[p++];
		}
	}
	delete[] paths;
}

void buildMST()
{
	for (int i = 0; i < n; i++)
	{
		key[i] = false;
		mstSet[i] = false;
		result[i] = -1;
	}
	mstSet[0] = true;
	result[0] = -1;
	for (int i = 1; i < n; i++)
	{
		runBuildKey(key, mstSet);
		runMinKey(key, mstSet);
		int u = 0;
		int ji = 0;
		float mmKey = FLT_MAX;
		for (int j = 0; j < procCount; j++) {
			if (mm[j] < mmKey) {
				u = minKeys[j];
				mmKey = mm[j];
				ji = jj[j];
			}
		}
		result[u] = ji;
		mstSet[u] = true;
	}
	delete[] key;
	delete[] mstSet;
}

void runBuildKey(bool* key, bool* mstSet)
{

	//int oStartIndex = startIndex;
	//int oCount = count;
	MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	//MPI_Send(&oStartIndex, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	//MPI_Send(&oCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	//MPI_Send(&key[oStartIndex], count, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
	//MPI_Send(&mstSet[oStartIndex], count, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
	//threads[i] = thread(buildKey, key, mstSet, startIndex, count);
	if (count != nominalCount)
	{
		break;
	}

	for (int i = 1; i < procCount; i++)
	{
		int t;
		MPI_Status status;
		MPI_Recv(&t, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
	}
}

void runMinKey(bool* key, bool* mstSet)
{
	int nominalCount = n / (procCount - 1);
	for (int startIndex = 0, i = 1; startIndex < n; startIndex += nominalCount, i++)
	{
		int count = nominalCount;
		if (startIndex + nominalCount < n && startIndex + nominalCount * 2 >= n)
		{
			count = n - startIndex;
		}
		int oper = OPER_MIN;
		int oStartIndex = startIndex;
		int oCount = count;
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&oStartIndex, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&oCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&key[oStartIndex], count, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
		MPI_Send(&mstSet[oStartIndex], count, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
		//threads[i] = thread(minKey, ref(minKeys[i]), ref(mm[i]), ref(jj[i]), key, mstSet, startIndex, count);
		if (count != nominalCount)
		{
			break;
		}
	}
	for (int i = 0; i < procCount; i++)
	{
		//threads[i].join();
	}
}

void buildKey(int pfrom, int count)
{
	int k = 0;
	for (int i = pfrom; k < count; i++, k++)
	{
		if (mstSet[i])
		{
			key[i] = false;
		}
		else
		{
			bool f = false;
			for (int j = 0; j < n; j++)
			{
				if (!mstSet[j])
				{
					continue;
				}
				if (graph[i * n + j] > 0)
				{
					f = true;
					break;
				}
			}
			key[i] = f;
		}
	}
}

minResp minKey(int pfrom, int count)
{
	float min = FLT_MAX;
	int minIndex = 0;
	int jetIndex = 0;
	int k = 0;
	for (int i = pfrom; k < count; i++, k++)
	{
		if (key[i])
		{
			for (int j = 0; j < n; j++)
			{
				if (mstSet[j])
				{
					float t = graph[i * n + j];
					if (t < min)
					{
						min = t;
						minIndex = i;
						jetIndex = j;
					}
				}
			}
		}
	}
	minResp result;
	result.v = minIndex;
	result.m = min;
	result.pj = jetIndex;
}

void writeMST()
{
	cout << "Writing mst ..." << endl;
	ofstream outf("mst.bin", ios::out | ios::binary);
	if (!outf)
	{
		cerr << "MST writing error!" << endl;
		exit(1);
	}
	outf.write((char*)(&n), sizeof(n));
	outf.write((char*)result, n * sizeof(int));
	outf.close();
}

void release()
{
	delete[] graph;
	delete[] result;
	delete[] points;
	delete[] minKeys;
	delete[] mm;
	delete[] jj;
	delete[] key;
	delete[] mstSet;
	cout << "Completed" << endl;
}