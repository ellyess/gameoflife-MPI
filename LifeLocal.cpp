#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include <time.h>
#include <fstream>

using namespace std;

// defining the game of life settings
const int rows = 100;
const int cols = 100;
const int period = 100;
bool periodic = true;

// initialising variables
int id, p, tag_num = 1;
int p_rows, p_cols;
int id_row, id_column;
int pad_rows, pad_cols;
int sub_rows, sub_cols;

// initialising variables for MPI Datatypes
vector<vector<int>> sub_length(16);
vector<vector<MPI_Aint>> addresses(16);
vector<vector<MPI_Datatype>> typelist(16);
vector<MPI_Datatype> temp_type(16);
vector<vector<bool>> recv_data(16);

// function that locates if a neighbour is alive
int neighbours(bool state[], int i, int j)
{
	int count = 0;
	// left
	if (state[i * pad_cols + j - 1] == true) count++;
	// right
	if (state[i * pad_cols + j + 1] == true) count++;
	// up
	if (state[(i - 1) * pad_cols + j] == true) count++;
	// down
	if (state[(i + 1) * pad_cols + j] == true) count++;
	// up left
	if (state[(i - 1) * pad_cols + j - 1] == true) count++;
	// up right
	if (state[(i - 1) * pad_cols + j + 1] == true) count++;
	// down left
	if (state[(i + 1) * pad_cols + j - 1] == true) count++;
	// down right
	if (state[(i + 1) * pad_cols + j + 1] == true) count++;
	return count;
}

// create a padded version of the targeted matrix
void pad(bool sub_grid[], bool pad_grid[], int sub_rows, int sub_cols)
{
	// copying the initial state into the larger padded array
	for (int i = 0; i < sub_rows; i++) {
		for (int j = 0; j < sub_cols; j++) {
			pad_grid[(i + 1) * (sub_cols + 2) + j + 1] = sub_grid[i * sub_cols + j];
		}
	}
}

// the core function to play game of life
void life(bool sub_grid[], bool pad_grid[])
{
	pad(sub_grid, pad_grid, sub_rows, sub_cols);

	// looping through the matrix and counting neighbour to determine fate of cell
	for (int i = 1; i < pad_rows - 1; i++) {
		for (int j = 1; j < pad_cols - 1; j++) {
			int total = 0;

			// conditions for an alive cell
			if (pad_grid[i * pad_cols + j] == true) {
				total = neighbours(pad_grid, i, j);
				if ((2 <= total) && (total <= 3))
					sub_grid[(i - 1) * sub_cols + (j - 1)] = true;
				else
					sub_grid[(i - 1) * sub_cols + (j - 1)] = false;
			}
			// conditions for a dead cell
			else {
				total = neighbours(pad_grid, i, j);
				if (total == 3)
					sub_grid[(i - 1) * sub_cols + (j - 1)] = true;
				else
					sub_grid[(i - 1) * sub_cols + (j - 1)] = false;
			}
		}
	}
}

// this function calculates the most effecient way to split the grids over
// the cores
void find_dimensions(int p, int& p_rows, int& p_cols, int& sub_rows, int& sub_cols)
{

	int min_gap = max(rows, cols);
	for (int i = 1; i <= p; i++)
		if (p % i == 0)
		{
		int gap = abs(cols / (p / i) - rows / i);
			if (gap < min_gap)
			{
				min_gap = gap;
				p_rows = i;
				p_cols = p / i;
			}
		}

	// calculating the number of rows for each processor subgrid
	int row_remainder = rows % p_rows;
	sub_rows = rows / p_rows;
	if (id < row_remainder * p_cols)
		sub_rows++;

	// calculating the number of columns for each processor subgrid
	int col_remainder = cols % p_cols;
	sub_cols = cols / p_cols;
	for (int i = 0; i < col_remainder; i++)
		if (id % p_cols == i)
			sub_cols++;

	if (id == 0)
		cout << "Divide " << p << " into " << p_rows << " by " << p_cols << " grid" << endl;

	cout << "processor " << id << " block is " << sub_rows << " x " << sub_cols << endl;
}

// indexing the processor id
void id_to_index(int id, int& id_row, int& id_column)
{
	id_column = id % p_cols;
	id_row = id / p_cols;
}

// finding the id from the index
int id_from_index(int id_row, int id_column, bool periodic)
{
	// if periodic it will find the id of of the processor opposite if at boundary
	if (periodic == true)
		return ((id_row + p_rows) % p_rows) * p_cols + ((id_column + p_cols) % p_cols);

	// else no processor exist next to boundary
	else {
		if (id_row >= p_rows || id_row < 0)
			return -1;
		if (id_column >= p_cols || id_column < 0)
			return -1;

		return id_row * p_cols + id_column;
	}
}

// this function creates data types for each processor to understand where to send and recieve
// data from their neighbours and also does the communication for periodic and non-periodic
// boundaries
void MPI_Communication(bool *pad_grid, int pad_cols, int pad_rows, int id)
{
	// top left send
	sub_length[0].push_back(1);
	MPI_Aint temp0;
	MPI_Get_address(&pad_grid[pad_cols + 1], &temp0);
	addresses[0].push_back(temp0);
	typelist[0].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[0].size(), &sub_length[0][0], &addresses[0][0], &typelist[0][0], &temp_type[0]);
	MPI_Type_commit(&temp_type[0]);

	// top middle send
	sub_length[1].push_back(pad_cols - 2);
	MPI_Aint temp1;
	MPI_Get_address(&pad_grid[pad_cols + 1], &temp1);
	addresses[1].push_back(temp1);
	typelist[1].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[1].size(), &sub_length[1][0], &addresses[1][0], &typelist[1][0], &temp_type[1]);
	MPI_Type_commit(&temp_type[1]);

	// top right send
	sub_length[2].push_back(1);
	MPI_Aint temp2;
	MPI_Get_address(&pad_grid[pad_cols + pad_cols - 2], &temp2);
	addresses[2].push_back(temp2);
	typelist[2].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[2].size(), &sub_length[2][0], &addresses[2][0], &typelist[2][0], &temp_type[2]);
	MPI_Type_commit(&temp_type[2]);

	// left send
	for (int i = 1; i < pad_rows - 1; i++)
	{
		sub_length[3].push_back(1);
		MPI_Aint temp3;
		MPI_Get_address(&pad_grid[i * pad_cols + 1], &temp3);
		addresses[3].push_back(temp3);
		typelist[3].push_back(MPI_C_BOOL);
	}
	MPI_Type_create_struct(sub_length[3].size(), &sub_length[3][0], &addresses[3][0], &typelist[3][0], &temp_type[3]);
	MPI_Type_commit(&temp_type[3]);

	// right send
	for (int i = 1; i < pad_rows - 1; i++)
	{
		sub_length[4].push_back(1);
		MPI_Aint temp4;
		MPI_Get_address(&pad_grid[i * pad_cols + (pad_cols - 2)], &temp4);
		addresses[4].push_back(temp4);
		typelist[4].push_back(MPI_C_BOOL);
	}
	MPI_Type_create_struct(sub_length[4].size(), &sub_length[4][0], &addresses[4][0], &typelist[4][0], &temp_type[4]);
	MPI_Type_commit(&temp_type[4]);

	// bottom left send
	sub_length[5].push_back(1);
	MPI_Aint temp5;
	MPI_Get_address(&pad_grid[((pad_rows - 2) * pad_cols) + 1], &temp5);
	addresses[5].push_back(temp5);
	typelist[5].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[5].size(), &sub_length[5][0], &addresses[5][0], &typelist[5][0], &temp_type[5]);
	MPI_Type_commit(&temp_type[5]);

	//  bottom middle send
	sub_length[6].push_back(pad_cols - 2);
	MPI_Aint temp6;
	MPI_Get_address(&pad_grid[((pad_rows - 2) * pad_cols) + 1], &temp6);
	addresses[6].push_back(temp6);
	typelist[6].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[6].size(), &sub_length[6][0], &addresses[6][0], &typelist[6][0], &temp_type[6]);
	MPI_Type_commit(&temp_type[6]);

	// bottom right send
	sub_length[7].push_back(1);
	MPI_Aint temp7;
	MPI_Get_address(&pad_grid[((pad_rows - 2) * pad_cols) + (pad_cols - 2)], &temp7);
	addresses[7].push_back(temp7);
	typelist[7].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[7].size(), &sub_length[7][0], &addresses[7][0], &typelist[7][0], &temp_type[7]);
	MPI_Type_commit(&temp_type[7]);

	// top left receive
	sub_length[8].push_back(1);
	MPI_Aint temp8;
	MPI_Get_address(&pad_grid[0], &temp8);
	addresses[8].push_back(temp8);
	typelist[8].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[8].size(), &sub_length[8][0], &addresses[8][0], &typelist[8][0], &temp_type[8]);
	MPI_Type_commit(&temp_type[8]);

	// top middle receive
	sub_length[9].push_back(pad_cols - 2);
	MPI_Aint temp9;
	MPI_Get_address(&pad_grid[1], &temp9);
	addresses[9].push_back(temp9);
	typelist[9].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[9].size(), &sub_length[9][0], &addresses[9][0], &typelist[9][0], &temp_type[9]);
	MPI_Type_commit(&temp_type[9]);

	// top right receive
	sub_length[10].push_back(1);
	MPI_Aint temp10;
	MPI_Get_address(&pad_grid[pad_cols - 1], &temp10);
	addresses[10].push_back(temp10);
	typelist[10].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[10].size(), &sub_length[10][0], &addresses[10][0], &typelist[10][0], &temp_type[10]);
	MPI_Type_commit(&temp_type[10]);

	// left receive
	for (int i = 1; i < pad_rows - 1; i++)
	{
		sub_length[11].push_back(1);
		MPI_Aint temp11;
		MPI_Get_address(&pad_grid[i * pad_cols], &temp11);
		addresses[11].push_back(temp11);
		typelist[11].push_back(MPI_C_BOOL);
	}
	MPI_Type_create_struct(sub_length[11].size(), &sub_length[11][0], &addresses[11][0], &typelist[11][0], &temp_type[11]);
	MPI_Type_commit(&temp_type[11]);

	// right receive
	for (int i = 1; i < pad_rows - 1; i++)
	{
		sub_length[12].push_back(1);
		MPI_Aint temp12;
		MPI_Get_address(&pad_grid[i * pad_cols + (pad_cols - 1)], &temp12);
		addresses[12].push_back(temp12);
		typelist[12].push_back(MPI_C_BOOL);
	}
	MPI_Type_create_struct(sub_length[12].size(), &sub_length[12][0], &addresses[12][0], &typelist[12][0], &temp_type[12]);
	MPI_Type_commit(&temp_type[12]);

	// bottom left receive
	sub_length[13].push_back(1);
	MPI_Aint temp13;
	MPI_Get_address(&pad_grid[(pad_rows - 1) * pad_cols], &temp13);
	addresses[13].push_back(temp13);
	typelist[13].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[13].size(), &sub_length[13][0], &addresses[13][0], &typelist[13][0], &temp_type[13]);
	MPI_Type_commit(&temp_type[13]);

	// bottom middle receive
	sub_length[14].push_back(pad_cols - 2);
	MPI_Aint temp14;
	MPI_Get_address(&pad_grid[((pad_rows - 1) * pad_cols) + 1], &temp14);
	addresses[14].push_back(temp14);
	typelist[14].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[14].size(), &sub_length[14][0], &addresses[14][0], &typelist[14][0], &temp_type[14]);
	MPI_Type_commit(&temp_type[14]);

	// bottom right receive
	sub_length[15].push_back(1);
	MPI_Aint temp15;
	MPI_Get_address(&pad_grid[(pad_rows * pad_cols) - 1], &temp15);
	addresses[15].push_back(temp15);
	typelist[15].push_back(MPI_C_BOOL);
	MPI_Type_create_struct(sub_length[15].size(), &sub_length[15][0], &addresses[15][0], &typelist[15][0], &temp_type[15]);
	MPI_Type_commit(&temp_type[15]);

	int cnt_type_send = 0;
	MPI_Request* request = new MPI_Request[8 * 2];
	int cnt = 0;
	tag_num = 0;
	int prime_fix = 0;
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			int com_i = id_row + i;
			int com_j = id_column + j;
			int com_id = id_from_index(com_i, com_j, periodic);

			// this if statement allows the compiler to deal with prime number of processors
			if (prime_fix != 4)
			{
				if (com_id >= 0 && com_id < p)
				{
					// sending the data from one type and recieving from the corresponding type on the other processor
					MPI_Isend(MPI_BOTTOM, 1, temp_type[cnt_type_send], com_id, tag_num, MPI_COMM_WORLD, &request[cnt * 2]);
					MPI_Irecv(MPI_BOTTOM, 1, temp_type[cnt_type_send + 8], com_id, (8 - tag_num), MPI_COMM_WORLD, &request[cnt * 2 + 1]);
					cnt++;
				}
				cnt_type_send++;
			}
			tag_num++;
			prime_fix++;
		}
	}
	MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);
	for (int i = 0; i < 16; i++)
		MPI_Type_free(&temp_type[i]);
}

int main(int argc, char* argv[])
{
	// initialising MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 10);

	// starting a timer
	double start = MPI_Wtime();

	// calling the functions to identify the processors dimensions and indexing
	find_dimensions(p, p_rows, p_cols, sub_rows, sub_cols);
	id_to_index(id, id_row, id_column);
	id_from_index(id_row, id_column, periodic);

	// outputting the settings of game of life for postprocessing
	if (id == 0) {
		fstream settings;
		settings.open("settings.txt", fstream::out);
		settings << rows << "\r\n";
		settings << cols << "\r\n";
		settings << period << "\r\n";
		settings << p_rows << "\r\n";
		settings << p_cols << "\r\n";
		settings.close();
	}

	// creating the dimensions of the padded array
	pad_rows = sub_rows + 2;
	pad_cols = sub_cols + 2;

	// initialising a random initial state for Life
	auto sub_grid = new bool[(sub_rows * sub_cols)]{ false };
	double perc_alive = 0.5;
	for (int i = 0; i < sub_rows; i++) {
		for (int j = 0; j < sub_cols; j++) {
			double num_rand = ((double)rand() / double(RAND_MAX));
			if (num_rand > 1 - perc_alive) { sub_grid[i * sub_cols + j] = true; }
		}
	}

	// opening a file that will have each processors subgrid at each generation
	fstream outfile;
	outfile.open("output_" + to_string(id) + ".txt", fstream::out);


	// outputs the dimensions per processor currently to their own file each
	// if time permits make this write to one file
	fstream dimensions;
	dimensions.open("dimensions_" + to_string(id) + ".txt", fstream::out);
	dimensions << sub_rows << " " << sub_cols;
	dimensions.close();

	// initial loop for a number of generations
	int gen = 0;
	auto pad_grid = new bool[(pad_cols * pad_rows)];
	while (gen <= period) {

		// creating a larger array for padded
		for (int i = 0; i < pad_rows * pad_cols; i++)
			pad_grid[i] = false;

		// copying the initial state into the larger padded array
		pad(sub_grid, pad_grid, sub_rows, sub_cols);

		// printing to a file for each processor that will contain
		for (int i = 0; i < sub_rows; i++) {
			for (int j = 0; j < sub_cols; j++) {
				outfile << sub_grid[i * sub_cols + j] << " ";
			}
			outfile << "\r\n";
		}

		outfile << "\r\n";

		// allowing each processor to communicate and play game of life
		MPI_Communication(pad_grid, pad_cols, pad_rows, id);
		life(sub_grid, pad_grid);

		gen++;
	}

	// ending the timer
	if (id == 0) {
		double end = MPI_Wtime();
		cout << "The process took " << end - start << " seconds to run." << endl;
	}

	outfile.close();
	delete[] pad_grid;
	delete[] sub_grid;
	MPI_Finalize();
}
