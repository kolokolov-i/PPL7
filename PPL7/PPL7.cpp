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

struct minResp
{
	int v;
	float m;
	int pj;
};

const int OPER_UPDATE = 1;
const int OPER_BUILD = 2;
const int OPER_MIN = 3;
const int OPER_EXIT = 4;
const int OPER_INIT = 5;
const int OPER_UKEY = 6;

void mReadGraph();
void mWriteMST();
void mInitWorkers();
void release();
void mBuildMST();
void mRunBuildKey();
void mRunMinKey();
void wBuildKey();
minResp wMinKey();
void worker();

int procCount = 3;
int curRank;
int n; // vertex count
int pathCount;
float* graph;
int* result;
int* minKeys;
int* jj;
float* mm;
bool* key;
bool* mstSet;
int* partsSize;
int* starts;
int partFrom;
int partSize;

int main(int argc, char* argv[])
{
	int ierr;
	ierr = MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) return ierr;
	MPI_Comm_rank(MPI_COMM_WORLD, &curRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	if (curRank == 0) {
		mReadGraph();
		minKeys = new int[procCount];
		mm = new float[procCount];
		jj = new int[procCount];
		key = new bool[n];
		mstSet = new bool[n];
		result = new int[n];
		mInitWorkers();
		mBuildMST();
		mWriteMST();
	}
	else {
		worker();
	}
	MPI_Finalize();
	release();
}

void mInitWorkers() {
	int nominalCount = n / (procCount - 1);
	partsSize = new int[procCount];
	starts = new int[procCount];
	for (int startIndex = 0, i = 1; startIndex < n; startIndex += nominalCount, i++) {
		int count = nominalCount;
		if (startIndex + nominalCount < n && startIndex + nominalCount * 2 >= n) {
			count = n - startIndex;
		}
		partsSize[i] = count;
		starts[i] = startIndex;
	}
	for (int i = 1; i < procCount; i++) {
		int oper = OPER_INIT;
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&starts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&partsSize[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(graph, n * n, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
	}
}

void mUpdateWorkers(int u) {
	int uu = u;
	for (int i = 1; i < procCount; i++) {
		int oper = OPER_UPDATE;
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&uu, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
}

void mUpdateKey() {
	for (int i = 1; i < procCount; i++) {
		int oper = OPER_UKEY;
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Status status;
		MPI_Recv(&key[starts[i]], partsSize[i], MPI_C_BOOL, i, 0, MPI_COMM_WORLD, &status);
	}
	for (int i = 1; i < procCount; i++) {
		MPI_Send(key, n, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
	}
}

void worker() {
	bool flag = true;
	while (flag) {
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
			for (int i = 0; i < n; i++)
			{
				key[i] = false;
				mstSet[i] = false;
			}
			mstSet[0] = true;
			break;
		case OPER_EXIT:
			return;
			break;
		case OPER_BUILD:
			wBuildKey();
			break;
		case OPER_MIN:
			minResp res = wMinKey();
			break;
		case OPER_UPDATE:
			int u;
			MPI_Recv(&u, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			mstSet[u] = true;
			break;
		case OPER_UKEY:
			MPI_Send(&key[partFrom], partSize, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(key, n, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, &status);
			break;
		}
	}
}

void mReadGraph()
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

void mBuildMST()
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
		mRunBuildKey();
		mUpdateKey();
		mRunMinKey();
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
		mUpdateWorkers(u);
	}
	delete[] key;
	delete[] mstSet;
}

void mRunBuildKey()
{
	int oper = OPER_BUILD;
	for (int i = 1; i < procCount; i++) {
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
}

void mRunMinKey()
{
	int oper = OPER_MIN;
	for (int i = 1; i < procCount; i++) {
		MPI_Send(&oper, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

	}
	for (int i = 1; i < procCount; i++) {

	}
}

void wBuildKey()
{
	int k = 0;
	for (int i = partFrom; k < partSize; i++, k++)
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

minResp wMinKey()
{
	float min = FLT_MAX;
	int minIndex = 0;
	int jetIndex = 0;
	int k = 0;
	for (int i = partFrom; k < partSize; i++, k++)
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
	return result;
}

void mWriteMST()
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
	delete[] minKeys;
	delete[] mm;
	delete[] jj;
	delete[] key;
	delete[] mstSet;
	cout << "Completed" << endl;
}