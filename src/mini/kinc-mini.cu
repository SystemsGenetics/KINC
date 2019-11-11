#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../cuda/linalg.cu"
#include "../cuda/sort.cu"
#include "../cuda/fetchpair.cu"
#include "../cuda/outlier.cu"
#include "../cuda/gmm.cu"
#include "../cuda/pearson.cu"
#include "../cuda/spearman.cu"
#include "../cuda/similarity.cu"



typedef enum ClusteringMethod ClusteringMethod;
typedef enum CorrelationMethod CorrelationMethod;



struct PairwiseIndex
{
	int x;
	int y;
};



PairwiseIndex pairwise_index(size_t index)
{
	// compute pairwise index from scalar index
	size_t pos {0};
	size_t x {0};

	while ( pos + x <= index )
	{
		pos += x;
		++x;
	}

	// return pairwise index
	return {
		static_cast<int>(x),
		static_cast<int>(index - pos)
	};
}



void operator++(PairwiseIndex& index)
{
	if ( ++index.y >= index.x )
	{
		index.y = 0;
		++index.x;
	}
}



struct Pair
{
	char K;
	std::vector<char> labels;
	std::vector<float> correlations;
};



#define CUDA_SAFE_CALL(ans) check((ans), #ans, __FILE__, __LINE__)

inline void check(cudaError_t err, const char *func, const char *file, const int line)
{
	if ( err != cudaSuccess ) {
		fprintf(stderr, "CUDA error at %s:%d\n", file, line);
		fprintf(stderr, "%s %s\n", cudaGetErrorString(err), func);
		exit(-1);
	}
}



template<typename T>
T * CUDABuffer(size_t size)
{
	T *ptr;
	CUDA_SAFE_CALL(cudaMallocManaged((void **)&ptr, size * sizeof(T)));

	return ptr;
}



void load_dataframe(const char *filename, int *p_rows, int *p_cols, float **p_data)
{
	// create input stream
	std::ifstream in(filename);

	// read sample names from first line
	std::string line;
	std::getline(in, line);

	// determine number of samples
	int cols {0};
	std::stringstream ss(line);

	while ( !ss.eof() )
	{
		std::string colname;
		ss >> colname;

		cols++;
	}

	// read data from input file
	int rows {0};
	std::vector<float> values;

	while ( !in.eof() )
	{
		// read a line from the input file
		std::getline(in, line);

		std::stringstream ss(line);

		// read row name
		std::string rowname;
		ss >> rowname;

		// read data elements
		for ( int i = 0; i < cols; i++ )
		{
			std::string token;
			ss >> token;

			// if token matches nan token then it as such
			if ( token == "NA" )
			{
				values.push_back(NAN);
			}

			// else this is a normal floating point expression
			else
			{
				float value;
				ss >> value;

				values.push_back(value);
			}
		}

		// increment number of rows
		rows++;
	}

	// initialize dataframe
	float *data = CUDABuffer<float>(rows * cols);

	memcpy(data, values.data(), rows * cols * sizeof(float));

	// save outputs
	*p_rows = rows;
	*p_cols = cols;
	*p_data = data;
}



int nextPowerTwo(int n)
{
	int pow2 = 2;
	while ( pow2 < n )
	{
		pow2 *= 2;
	}

	return pow2;
}



template<class T>
std::vector<T> makeVector(const T* data, int size)
{
	std::vector<T> v(size);

	memcpy(v.data(), data, size * sizeof(T));
	return v;
}



int main(int argc, char **argv)
{
	// parse command-line arguments
	if ( argc != 2 )
	{
		fprintf(stderr, "usage: ./kinc-mini <infile>\n");
		exit(-1);
	}

	const char *filename = argv[1];

	// load dataframe
	int geneSize;
	int sampleSize;
	float *expressions;

	load_dataframe(filename, &geneSize, &sampleSize, &expressions);

	printf("loaded dataframe (%d x %d)\n", geneSize, sampleSize);

	// copy dataframe to GPU
	CUDA_SAFE_CALL(cudaMemPrefetchAsync(expressions, geneSize * sampleSize, 0, 0));

	// initialize execution parameters
	int globalWorkSize {4096};
	int localWorkSize {32};
	int minSamples {30};
	float minExpression {-INFINITY};
	ClusteringMethod clusMethod {ClusteringMethod_GMM};
	CorrelationMethod corrMethod {CorrelationMethod_Spearman};
	int minClusters {1};
	int maxClusters {5};
	Criterion criterion {ICL};
	bool removePreOutliers {true};
	bool removePostOutliers {true};

	// initialize buffers
	int W {globalWorkSize};
	int N {sampleSize};
	int N_pow2 {nextPowerTwo(N)};
	int K {maxClusters};

	int2 *	in_index			= CUDABuffer<int2>(1 * W);
	float *  work_x			  = CUDABuffer<float>(N_pow2 * W);
	float *  work_y			  = CUDABuffer<float>(N_pow2 * W);
	float2 * work_gmm_data		= CUDABuffer<float2>(N * W);
	char *	work_gmm_labels	 = CUDABuffer<char>(N * W);
	float *  work_gmm_pi		 = CUDABuffer<float>(K * W);
	float2 * work_gmm_mu		 = CUDABuffer<float2>(K * W);
	float4 * work_gmm_sigma	  = CUDABuffer<float4>(K * W);
	float4 * work_gmm_sigmaInv	= CUDABuffer<float4>(K * W);
	float *  work_gmm_normalizer = CUDABuffer<float>(K * W);
	float2 * work_gmm_MP		 = CUDABuffer<float2>(K * W);
	int *	work_gmm_counts	 = CUDABuffer<int>(K * W);
	float *  work_gmm_logpi	  = CUDABuffer<float>(K * W);
	float *  work_gmm_gamma	  = CUDABuffer<float>(N * K * W);
	char *	out_K				= CUDABuffer<char>(1 * W);
	char *	out_labels		  = CUDABuffer<char>(N * W);
	float *  out_correlations	= CUDABuffer<float>(K * W);

	// initialize output
	std::vector<Pair> pairs;

	// iterate through all pairs
	int workBlockStart {0};
	int workBlockSize {32768};
	PairwiseIndex index = pairwise_index(workBlockStart);

	for ( int i = 0; i < workBlockSize; i += globalWorkSize )
	{
		printf("%8d %4d %4d\n", workBlockStart, index.x, index.y);

		// write input buffers to device
		int numPairs {min(globalWorkSize, workBlockSize - i)};

		for ( int j = 0; j < numPairs; ++j )
		{
			in_index[j] = { index.x, index.y };
			++index;
		}

		CUDA_SAFE_CALL(cudaMemPrefetchAsync(in_index, W * sizeof(int2), 0, 0));

		// execute similiarity kernel
		Similarity_compute<<<globalWorkSize, localWorkSize>>>(
			clusMethod,
			corrMethod,
			removePreOutliers,
			removePostOutliers,
			numPairs,
			expressions,
			sampleSize,
			in_index,
			minExpression,
			minSamples,
			minClusters,
			maxClusters,
			criterion,
			work_x,
			work_y,
			work_gmm_data,
			work_gmm_labels,
			work_gmm_pi,
			work_gmm_mu,
			work_gmm_sigma,
			work_gmm_sigmaInv,
			work_gmm_normalizer,
			work_gmm_MP,
			work_gmm_counts,
			work_gmm_logpi,
			work_gmm_gamma,
			out_K,
			out_labels,
			out_correlations
		);
		CUDA_SAFE_CALL(cudaGetLastError());

		// read results from device
		CUDA_SAFE_CALL(cudaMemPrefetchAsync(out_K, W * sizeof(char), cudaCpuDeviceId, 0));
		CUDA_SAFE_CALL(cudaMemPrefetchAsync(out_labels, W * N * sizeof(char), cudaCpuDeviceId, 0));
		CUDA_SAFE_CALL(cudaMemPrefetchAsync(out_correlations, W * K * sizeof(float), cudaCpuDeviceId, 0));

		// wait for everything to finish
		CUDA_SAFE_CALL(cudaStreamSynchronize(0));

		// save results
		for ( int j = 0; j < numPairs; ++j )
		{
			// get pointers to the cluster labels and correlations for this pair
			const char *labels = &out_labels[j * sampleSize];
			const float *correlations = &out_correlations[j * maxClusters];

			Pair pair;

			// save the number of clusters
			pair.K = out_K[j];

			// save the cluster labels and correlations (if the pair was able to be processed)
			if ( pair.K > 0 )
			{
				pair.labels = makeVector(labels, sampleSize);
				pair.correlations = makeVector(correlations, maxClusters);
			}

			pairs.push_back(pair);
		}
	}

	return 0;
}
