#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
/*
	Current version of the algorithm takes 2 strings not longer than 40 characters.
	Gap and mismatch penalties are 0
	Match score is 1
	Algorithm returns the longest substring by tracing back the score matrix that scored with the penalties/scores above
	The match/mismatch/gap values can be adjusted to get the prefered result
*/

/*	
	TODO
	1) Handle the exception where the size is smaller than shorter string's lenght
		Distribute more cells to process to each thread
		Handle the offsets and number of data sent to each thread
	2) Clean up comments
	3) Add match/mismatch/gap values to the code instead of using magic numbers
	[Optional] Check out if there is a better way of broadcasting
	[Optional] Check out if there is a better way of scattering
*/
using namespace std;

class SmithWaterman {
private:

	#pragma region Variables
	char str1[40], str2[40];
	int	s1Lenght, s2Lenght;
	int diagonalCount, vectorRowCount, minCount, maxCount;
	vector<vector<uint8_t>> scoreMatrix;
	uint8_t* currentCells;
	vector<int> cellsPerRank;
	vector<int> offsetCase1;
	vector<int> offsetCase2;
	vector<int> dataPerProcess;
	uint8_t localCellScore;
	uint8_t localRequired[3];
	#pragma endregion	

public:

	#pragma region Class Methods

	int getDiagCount() {

		return vectorRowCount;

	}

	SmithWaterman(string s1,
		string s2,
		int size) {

		strcpy(str1, s1.c_str());
		strcpy(str2, s2.c_str());
		s1Lenght = s1.size();
		s2Lenght = s2.size();

		diagonalCount = s1Lenght + s2Lenght - 1;
		vectorRowCount = diagonalCount + 2;
		minCount = min(s1Lenght, s2Lenght);
		maxCount = max(s1Lenght, s2Lenght);

		scoreMatrix.resize(vectorRowCount);
		offsetCase1.resize(size);
		offsetCase2.resize(size);
		dataPerProcess.resize(size, 0);

		for (int i = 0; i < size; i++) {
			offsetCase1[i] = i + 1;
			offsetCase2[i] = i;
		}

		for (int i = 0; i < vectorRowCount; i++) {
			scoreMatrix[i] = vector<uint8_t>(currentSize(i), 0);
		}

	}

	#pragma endregion	

	#pragma region MPI Calls
	void broadcast() {

		MPI_Bcast(&s1Lenght, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&s2Lenght, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&diagonalCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&vectorRowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&minCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(str1, s1Lenght, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(str2, s2Lenght, MPI_CHAR, 0, MPI_COMM_WORLD);

	}

	void scatter(int currentDiag,
		int size,
		int rank) {

		if (currentDiag <= (minCount + 1)) {
			dataPerProcess[currentDiag - 2] = 1;
		}
		else if (currentDiag > (maxCount + 1)) {
			dataPerProcess[currentSize(currentDiag)] = 0;
		}

		currentCells = scoreMatrix[currentDiag].data();

		//scattering local required and current diagonal 
		MPI_Scatterv(scoreMatrix[currentDiag - 2].data(), &dataPerProcess[0], currentDiag <= minCount + 1 ? &offsetCase2[0] : &offsetCase1[0], MPI_UINT8_T, &localRequired[0], dataPerProcess[rank], MPI_UINT8_T, 0, MPI_COMM_WORLD);	//diag
		MPI_Scatterv(scoreMatrix[currentDiag - 1].data(), &dataPerProcess[0], &offsetCase2[0], MPI_UINT8_T, &localRequired[1], dataPerProcess[rank], MPI_UINT8_T, 0, MPI_COMM_WORLD);	//left
		MPI_Scatterv(scoreMatrix[currentDiag - 1].data(), &dataPerProcess[0], &offsetCase1[0], MPI_UINT8_T, &localRequired[2], dataPerProcess[rank], MPI_UINT8_T, 0, MPI_COMM_WORLD);	//up


		//needed local variables
		//localRequired[0] = scoreMatrix[currentDiag - 2][currentDiag > minCount ? rank + 1 : rank];	//diagonal
		//localRequired[1] = scoreMatrix[currentDiag - 1][rank];										//left
		//localRequired[2] = scoreMatrix[currentDiag - 1][rank + 1];									//up
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void gather(int currentDiag, int rank) {

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(&localCellScore, dataPerProcess[rank], MPI_UINT8_T, &currentCells[0], &dataPerProcess[0], currentDiag <= minCount ? &offsetCase1[0] : &offsetCase2[0], MPI_UINT8_T, 0, MPI_COMM_WORLD);

	}

	#pragma endregion	

	#pragma region Algorithm Methods

	void score(int currentDiag,
		//uint8_t jobCoef,
		int rank) {

		//score for only active ranks
		//cout << "cD: " << currentDiag << " rank: " << rank << " local0: " << +localRequired[0] << " local1: " << +localRequired[1] << " local2: " << +localRequired[2] << endl;
		if (currentDiag <= minCount + 1) {
			localCellScore = str1[currentDiag - 2 - rank] == str2[rank] ? localRequired[0] + 1 : max(localRequired[1], localRequired[2]);
		}
		else {
			// current diagonal is smaller than min count
			localCellScore = str1[minCount - 1 - rank] == str2[currentDiag - minCount - 1 + rank] ? localRequired[0] + 1 : max(localRequired[1], localRequired[2]);
		}

	}

	void traceback() {

		uint8_t maxScore = 0;
		string result = "";
		int x, y;

		for (int i = 0; i < vectorRowCount; i++) {
			int col = currentSize(i);

			for (int j = 0; j < col; j++) {

				if (scoreMatrix[i][j] > maxScore) {
					maxScore = scoreMatrix[i][j];
					x = i;
					y = j;
				}

			}

		}

		while (scoreMatrix[x][y] != 0) {

			uint8_t left, up, diag, maxS;

			if (x <= minCount) {
				diag = scoreMatrix[x - 2][y - 1];
				left = scoreMatrix[x - 1][y - 1];
				up = scoreMatrix[x - 1][y];
			}
			else if (x == minCount + 1) {
				diag = scoreMatrix[x - 2][y];
				left = scoreMatrix[x - 1][y];
				up = scoreMatrix[x - 1][y + 1];
			}
			else {
				diag = scoreMatrix[x - 2][y + 1];
				left = scoreMatrix[x - 1][y];
				up = scoreMatrix[x - 1][y + 1];
			}

			maxS = max(diag, max(left, up));

			if (maxS == diag) {

				result = str1[x <= minCount ? y : (minCount - y - 1)] + result;

				if (x <= minCount + 1) {
					y = y - 1;
				}
				else {
					y = y + 1;
				}

				x = x - 2;
			}
			else if (maxS == left) {

				if (x <= minCount + 1) {
					y = y - 1;
				}

				x = x - 1;
			}
			else {

				if (x <= minCount + 1) {
					y = y;
				}
				else {
					y = y + 1;
				}

				x = x - 1;
			}

		}
		
		cout << result << endl;

	}

	#pragma endregion	

	#pragma region Helper Methods

	int currentSize(int currentDiag) {

		if (currentDiag <= minCount) return currentDiag + 1;
		if (currentDiag <= maxCount) return minCount + 1;

		return (minCount + 1) - (currentDiag - (maxCount));

	}

	int processCount(int currentDiag) {

		if (currentDiag <= minCount) return currentDiag - 1;
		if (currentDiag <= maxCount) return minCount;

		return (minCount + 1) - (currentDiag - (maxCount));

	}

	#pragma endregion	

};

int main(int argc, char* argv[]) {

	int rank, size;
	string s1 = "TAGTCACG";
	string s2 = "AGACTGTC";
	string s1f, s2f;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*
	if shorter string's length is larger than size, we need to distribute cells

	int maxProcessAtaTime = max(s1.size(), s2.size());
	int len_size = ceil(maxProcessAtaTime / size);
	int* lenght = new int[len_size];
	*/
	
	if (s1.length() <= s2.length()) {
		s1f = s1;
		s2f = s2;
	} else {
		s1f = s2;
		s2f = s1;
	}

	SmithWaterman Test(s1f, s2f, size);
	Test.broadcast();
	int diagCount = Test.getDiagCount();

	for (int i = 2; i < diagCount; i++) {

		MPI_Barrier(MPI_COMM_WORLD);
		Test.scatter(i, size, rank);

		if (rank < Test.processCount(i)) {
			Test.score(i, rank);
		}
			
		Test.gather(i,rank);

	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		Test.traceback();
	}		

	MPI_Finalize();

}
