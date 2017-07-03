#include "utility.h"
// /**
//         print error message then exit program
// */
// void outError(const char *error){
//         cerr << "ERROR: " << error << endl;
//         exit(EXIT_FAILURE);
// }
// 
// /**
//         print error message then exit program
// */
// void outError(string error){
//         outError(error.c_str());
// }
// 
// /**
//         print double error messages then exit program
// */
// void outError(const char *error, const char *msg){
//         string str = error;
//         str += msg;
//         outError(str);
// }
// 
// /**
//         print double error messages then exit program
// */
// void outError(const char *error, string msg){
//         string str = error;
//         str += msg;
//         outError(str);
// }
// 
// /**
//         convert string to integer, with error checking
//         @param str original string
//         @return the integer value
// */
// int convert_int(const char *str) throw (string){
//         char *endptr;
//         long i = strtol(str, &endptr, 10);
// 
//         if ((i == 0 && ((long) endptr - (const long) str) == 0) || abs(i) == HUGE_VALL || *endptr != 0) {
//                 string err = "Expecting integer, but found \"";
//                 err += str;
//                 err += "\" instead";
//                 throw err;
//         }
// 
//         return i;
// }/**
//         convert string to double, with error checking
//         @param str original string
//         @return the double
// */
// double convert_double(const char *str) throw (string){
//         char *endptr;
//         double d = strtod(str, &endptr);
//         if ((d == 0.0 && ((long) endptr - (const long) str) == 0) || fabs(d) == HUGE_VALF || *endptr != 0) {
//                 string err = "Expecting floating-point number, but found \"";
//                 err += str;
//                 err += "\" instead";
//                 throw err;
//         }
//         return d;
// }
// 
// string convertIntToString(int number)
// {
//    stringstream ss;//create a stringstream
//    ss << number;//add number to the stream
//    return ss.str();//return a string with the contents of the stream
// }
// 
// // string convertToLower(const string str)
// // {
// // 	string temp(str);
// // 	boost::algorithm::to_lower(temp);
// // 	return temp;
// // 	
// // }
// 
// int searchString(const string str, const StringVec& str_vector)
// {
// 	int size = str_vector.size();
// 	for ( int i = 0; i < size; i++ )
// 		if ( str == str_vector[i] ) return i;
// 	return -1;
// }
// int ordinaryRounding(double a)
// {
// 	return floor(a+0.5);
// }
// 

/* helpers for memory management */
double** allocDoubleMatrix(int rows, int cols) {
    double** matrix = (double**) calloc(rows, sizeof(double*));
    int i;
    for (i=0; i<rows; i++) {
        matrix[i] = (double*) calloc(cols, sizeof(double));
    }
    return(matrix);
}

void freeDoubleMatrix(double** matrix, int rows) {
    int i;
    for (i=0; i<rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// 
// int** allocIntMatrix(int rows, int cols) {
//     int** matrix = (int**) calloc(rows, sizeof(int*));
//     int i;
//     for (i=0; i<rows; i++) {
//         matrix[i] = (int*) calloc(cols, sizeof(int));
//     }
//     return(matrix);
// }
// 
// void freeIntMatrix(int** matrix, int rows) {
//     int i;
//     for (i=0; i<rows; i++) {
//         free(matrix[i]);
//     }
//     free(matrix);
// }
// 
 bool** allocBoolMatrix(int rows, int cols) {
     bool** matrix = (bool**) calloc(rows, sizeof(bool*));
     int i;
     for (i=0; i<rows; i++) {
         matrix[i] = (bool*) calloc(cols, sizeof(bool));
     }
     return(matrix);
 }
// 
 void freeBoolMatrix(bool** matrix, int rows) {
     int i;
     for (i=0; i<rows; i++) {
         free(matrix[i]);
     }
     free(matrix);
 }
// 
 double*** alloc3Ddouble(int dim1, int dim2, int dim3)
 {
     int i,j,k;
     double *** array = (double ***)malloc(dim1*sizeof(double**));
     
     for (i = 0; i< dim1; i++) {
         
         array[i] = (double **) malloc(dim2*sizeof(double *));
         
         for (j = 0; j < dim2; j++) {
             
             array[i][j] = (double *)malloc(dim3*sizeof(double));
         }
         
     }
     return array;
 }
// 
 void free3Ddouble(double*** array, int dim1, int dim2)
 {
     for ( int i=0; i < dim1; i++)
         freeDoubleMatrix(array[i], dim2);
     free(array);
 }
// 
// 
int argMax(double *a, const int N){
	double maximum=a[0];
	int argmax=0;
	for(int i=0;i<N;i++){
		if(maximum<a[i]){
			argmax=i;
			maximum=a[i];
		};
	};
	return argmax;
};
int argIntMax(int *a, const int N){
	int maximum=a[0];
	int argmax=0;
	for(int i=0;i<N;i++){
		if(maximum<a[i]){
			argmax=i;
			maximum=a[i];
		};
	};
	return argmax;
};
double Max(double *a, int N){
    double maximum=a[0];
    for(int i=0;i<N;i++){
        if(maximum<a[i]) maximum=a[i];
    };
    return maximum;
};
int intMax(int *a, int N){
    int maximum=a[0];
    for(int i=0;i<N;i++){
        if(maximum<a[i]) maximum=a[i];
    };
    return maximum;
};
double MaxMatrix(double **a, int N, int M){
    double maximum=a[0][0];
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            if(maximum<a[i][j]) maximum=a[i][j];
        };
    };
    return maximum;
};
int MaxIntMatrix(int** a, int N, int M)
{
    int maximum=a[0][0];
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            if(maximum<a[i][j]) maximum=a[i][j];
        };
    };
    return maximum;
}

double MaxDoubleMatrix(double** a, int N, int M)
{
    double maximum=a[0][0];
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            if(maximum<a[i][j]) maximum=a[i][j];
        };
    };
    return maximum;
}
// 
// /**
//  To read parameters for the density function, filename is the output of uniHMM for one chromosome (and at the moment one mark)
//  Format of the input file:
//  iteration    Pro.0   Pro.1   A.00    A.01    A.10    A.11    r0  r1  p0  p1  w0  w1
//  In the case: filename contains output for one mark, the density parameters for this mark will be store in densparam[0][]
//  */
// void readUniDensityParam(const char *filename, double **densparam)
// {
//     //start reading the data file
//     cout << "Reading model (density function) parameters from file: " << filename << ".........." << endl;
//     ostringstream err_str;
//     ifstream in;
//     int line_num = 0;
// 	try{
// 		// set the failbit and badbit
// 		in.exceptions(ios::failbit | ios::badbit);
//         // 			in.open(    markDataTable);
// 		in.open(filename);
// 		string line;
// 		string temp;
// 		// remove the failbit
// 		in.exceptions(ios::badbit);
// 		//read the first line: iteration    Pro.0   Pro.1   A.00    A.01    A.10    A.11    r0  r1  p0  p1  w0  w1
// 		getline(in, line);
// 		if (line == "")
// 		{
// 			err_str << "\tThe first line is empty!";
// 			throw err_str.str();
// 		}
// 		istringstream first_line(line);
//         int nCol = 13;
//         int toIgnore = nCol - 6;
//         for (int i = 0; i < nCol; i++)
//         {
//             if ( !(first_line >> temp) )
//             {
//                 err_str << "\tThe first line must contain at least " << nCol << " column_names";
//                 throw err_str.str();
//             }
//         }
//         for (; !in.eof(); line_num++)
//         {
// 			getline (in, line);
// 			if ( line == "" ) continue;
// 			istringstream line_in(line);
//             for (int i = 0; i < toIgnore; i++) //ignore the first 7 numbers
//             {
//                 if ( !(line_in >> temp) )
//                 {
//                     err_str << "Line " << line_num+2 << "contains less than " << nCol << " elements";
//                     throw err_str.str();
//                 }
//             }
//             for (int i = 0; i < 6; i++) //start reading the parameters of the density function
//             {
//                 if ( !(line_in >> temp) )
//                 {
//                     err_str << "Line " << line_num+2 << "contains less than " << nCol << " elements";
//                     throw err_str.str();
//                 }
//                 densparam[line_num][i] = convert_double(temp.c_str());
//             }
//         }
//         in.clear();
//         // set the failbit again
//         in.exceptions(ios::failbit | ios::badbit);
//         in.close();
//     }
//     catch (ios::failure)
// 	{
// 		outError("Could not read file: ", filename);
// 	}
//     cout << "Finish reading file" << filename << endl;
// }
// 
// /**
//     To read calls (0/1) from the uniHMM. At the moment input file is the output of uniHMM for one mark and one chromosome
//     Format of the input file:
//     Pr(state=0) Pr(state=1) maxState
//     (the number of lines of this file is equal to the number of bins of the chromosome plus one for the header line)
//     @param markID index of the histone mark. enrich is of dimension T*N where T is the number of bins, N is the number of marks. So markID indicates which column (or mark) in
//     enrich should we fetch the calls in.
//  */
// void readUniCall(const char *filename, bool **enrich, const int markID)
// {
//     //start reading the data file
//     cout << "Reading univariate calls from file: " << filename << "........" << endl;
//     ostringstream err_str;
//     ifstream in;
//     int line_num = 0;
// 	try{
// 		// set the failbit and badbit
// 		in.exceptions(ios::failbit | ios::badbit);
//         // 			in.open(params.markDataTable);
// 		in.open(filename);
// 		string line;
// 		string temp;
// 		// remove the failbit
// 		in.exceptions(ios::badbit);
// 		//read the first line: Pr(state=0) Pr(state=1) maxState
// 		getline(in, line);
// 		if (line == "")
// 		{
// 			err_str << "\tThe first line is empty!";
// 			throw err_str.str();
// 		}
// 		istringstream first_line(line);
//         int nCol = 3;
//         int toIgnore = nCol - 1;
//         for (int i = 0; i < nCol; i++)
//         {
//             if ( !(first_line >> temp) )
//             {
//                 err_str << "\tThe first line must contain at least " << nCol << " column_names";
//                 throw err_str.str();
//             }
//         }
//         for (; !in.eof(); line_num++)
//         {
// 			getline (in, line);
// 			if ( line == "" ) continue;
// 			istringstream line_in(line);
//             for (int i = 0; i < toIgnore; i++) //ignore the first 2 numbers
//             {
//                 if ( !(line_in >> temp) )
//                 {
//                     err_str << "Line " << line_num+2 << "contains less than " << nCol << " elements";
//                     throw err_str.str();
//                 }
//             }
//             //read the last number (index of the state having the maximum posterior
//             if ( !(line_in >> temp) )
//             {
//                 err_str << "Line " << line_num+2 << "contains less than " << nCol << " elements";
//                 throw err_str.str();
//             }
// if (line_num == 10000) cout<<"error here"<<endl;
//             enrich[line_num][markID] = (convert_int(temp.c_str()) == 0)? false : true;
//         }
//         in.clear();
//         // set the failbit again
//         in.exceptions(ios::failbit | ios::badbit);
//         in.close();
//     }
//     catch (ios::failure)
// 	{
// 		outError("Could not read file: ", filename);
// 	}
//     cout << "Finish reading file" << filename << endl;
// }
// 
// 

/* log likelihood */

/* x matrix [T x N] */
/*double sumP(int begin, int end, double** x, double* norm, int N, int T){
    double s, sum;
    int q, i;
    s = 0.0;
    for (q=begin; q<end; q++) {
        sum = 0.0;
        for (i=0; i<N; i++) {
            sum += x[q][i] * norm[q];
        }
        s += log(1.0 / sum);
    }
    return s;
}*/
