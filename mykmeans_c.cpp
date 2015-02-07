#include "stdlib.h"
#include "mex.h"
#include "cv.h"  
#include "highgui.h"

#include <list>
#include <vector>

using namespace std;
using namespace cv;




#ifndef USE_SSE2
volatile bool USE_SSE2 = false;
#endif

static inline float mydistance(const float* a, const float* b, int n)
{
    int j = 0; float d = 0.f;
#if CV_SSE
    if( USE_SSE2 )
    {
        float CV_DECL_ALIGNED(16) buf[4];
        __m128 d0 = _mm_setzero_ps(), d1 = _mm_setzero_ps();

        for( ; j <= n - 8; j += 8 )
        {
            __m128 t0 = _mm_sub_ps(_mm_loadu_ps(a + j), _mm_loadu_ps(b + j));
            __m128 t1 = _mm_sub_ps(_mm_loadu_ps(a + j + 4), _mm_loadu_ps(b + j + 4));
            d0 = _mm_add_ps(d0, _mm_mul_ps(t0, t0));
            d1 = _mm_add_ps(d1, _mm_mul_ps(t1, t1));
        }
        _mm_store_ps(buf, _mm_add_ps(d0, d1));
        d = buf[0] + buf[1] + buf[2] + buf[3];
    }
    else
#endif
    {
        for( ; j <= n - 4; j += 4 )
        {
            float t0 = a[j] - b[j], t1 = a[j+1] - b[j+1], t2 = a[j+2] - b[j+2], t3 = a[j+3] - b[j+3];
            //d += t0*t0*w[j] + t1*t1*w[j+1]  + t2*t2*w[j+2]  + t3*t3*w[j+3] ;
			//d += fabs(t0) + fabs(t1) + fabs(t2) + fabs(t3);
			d += t0*t0 + t1*t1  + t2*t2  + t3*t3 ;
        }
    }

    for( ; j < n; j++ )
    {
        float t = a[j] - b[j];
        //d += t*t*w[j] ;
		//d += fabs(t);
		d+=t*t;
    }
    return d;
}




static void generateCentersPP(const cv::Mat& _data, cv::Mat& _out_centers,int K, cv::RNG& rng, int trials)
{
    int i, j, k, dims = _data.cols, N = _data.rows;
    const float* data = _data.ptr<float>(0);
    size_t step = _data.step/sizeof(data[0]);
    vector<int> _centers(K);
    int* centers = &_centers[0];
    vector<float> _dist(N*3);
    float* dist = &_dist[0], *tdist = dist + N, *tdist2 = tdist + N;
    double sum0 = 0;

    centers[0] = (unsigned)rng % N;

    for( i = 0; i < N; i++ )
    {
        dist[i] = mydistance(data + step*i, data + step*centers[0], dims);
        sum0 += dist[i];
    }
    
    for( k = 1; k < K; k++ )
    {
        double bestSum = DBL_MAX;
        int bestCenter = -1;

        for( j = 0; j < trials; j++ )
        {
            double p = (double)rng*sum0, s = 0;
            for( i = 0; i < N-1; i++ )
                if( (p -= dist[i]) <= 0 )
                    break;
            int ci = i;
            for( i = 0; i < N; i++ )
            {
                tdist2[i] = std::min(mydistance(data + step*i, data + step*ci,dims), dist[i]);
                s += tdist2[i];
            }
            
            if( s < bestSum )
            {
                bestSum = s;
                bestCenter = ci;
                std::swap(tdist, tdist2);
            }
        }
        centers[k] = bestCenter;
        sum0 = bestSum;
        std::swap(dist, tdist);
    }

    for( k = 0; k < K; k++ )
    {
        const float* src = data + step*centers[k];
        float* dst = _out_centers.ptr<float>(k);
        for( j = 0; j < dims; j++ )
            dst[j] = src[j];
    }
}


static void generateRandomCenter(const vector<Vec2f>& box, float* center, RNG& rng)
{
    size_t j, dims = box.size();
    float margin = 1.f/dims;
    for( j = 0; j < dims; j++ )
        center[j] = ((float)rng*(1.f+margin*2.f)-margin)*(box[j][1] - box[j][0]) + box[j][0];
}






double mykmeans( cv::InputArray _data, int K, cv::InputOutputArray _bestLabels,cv::TermCriteria criteria, const float *w, int attempts, int flags, OutputArray _centers ){
    const int SPP_TRIALS = 3;
    cv::Mat data = _data.getMat();
    int N = data.rows > 1 ? data.rows : data.cols;
    int dims = (data.rows > 1 ? data.cols : 1)*data.channels();
    int type = data.depth();


    attempts = std::max(attempts, 1);
    CV_Assert( data.dims <= 2 && type == CV_32F && K > 0 );
	CV_Assert( data.dims <= 2  && K > 0 );
    _bestLabels.create(N, 1, CV_32S, -1, true);
    
    cv::Mat _labels, best_labels = _bestLabels.getMat();
    if( flags & CV_KMEANS_USE_INITIAL_LABELS )
    {
        CV_Assert( (best_labels.cols == 1 || best_labels.rows == 1) &&
                  best_labels.cols*best_labels.rows == N &&
                  best_labels.type() == CV_32S &&
                  best_labels.isContinuous());
        best_labels.copyTo(_labels);
    }
    else
    {
        if( !((best_labels.cols == 1 || best_labels.rows == 1) &&
             best_labels.cols*best_labels.rows == N &&
            best_labels.type() == CV_32S &&
            best_labels.isContinuous()))
            best_labels.create(N, 1, CV_32S);
        _labels.create(best_labels.size(), best_labels.type());
    }
    int* labels = _labels.ptr<int>();

    cv::Mat centers(K, dims, type), old_centers(K, dims, type);
    vector<int> counters(K);
    vector<Vec2f> _box(dims);
    Vec2f* box = &_box[0];

    double best_compactness = DBL_MAX, compactness = 0;
    RNG& rng = theRNG();
    int a, iter, i, j, k;

    if( criteria.type & TermCriteria::EPS )
        criteria.epsilon = std::max(criteria.epsilon, 0.);
    else
        criteria.epsilon = FLT_EPSILON;
    criteria.epsilon *= criteria.epsilon;

    if( criteria.type & TermCriteria::COUNT )
        criteria.maxCount = std::min(std::max(criteria.maxCount, 2), 100);
    else
        criteria.maxCount = 100;

    if( K == 1 )
    {
        attempts = 1;
        criteria.maxCount = 2;
    }

    const float* sample = data.ptr<float>(0);
    for( j = 0; j < dims; j++ )
        box[j] = Vec2f(sample[j], sample[j]);

    for( i = 1; i < N; i++ )
    {
        sample = data.ptr<float>(i);
        for( j = 0; j < dims; j++ )
        {
            float v = sample[j];
            box[j][0] = std::min(box[j][0], v);
            box[j][1] = std::max(box[j][1], v);
        }
    }

    for( a = 0; a < attempts; a++ )
    {
        double max_center_shift = DBL_MAX;
        for( iter = 0; iter < criteria.maxCount && max_center_shift > criteria.epsilon; iter++ )
        {
            swap(centers, old_centers);

            if( iter == 0 && (a > 0 || !(flags & KMEANS_USE_INITIAL_LABELS)) )
            {
                if( flags & KMEANS_PP_CENTERS )
                    generateCentersPP(data, centers, K, rng, SPP_TRIALS);
                else
                {
                    for( k = 0; k < K; k++ )
                        generateRandomCenter(_box, centers.ptr<float>(k), rng);
                }
            }
            else
            {
                if( iter == 0 && a == 0 && (flags & KMEANS_USE_INITIAL_LABELS) )
                {
                    for( i = 0; i < N; i++ )
                        CV_Assert( (unsigned)labels[i] < (unsigned)K );
                }
            
                // compute centers
                centers = Scalar(0);
                for( k = 0; k < K; k++ )
                    counters[k] = 0;

                for( i = 0; i < N; i++ )
                {
                    sample = data.ptr<float>(i);
                    k = labels[i];
                    float* center = centers.ptr<float>(k);
                    for( j = 0; j <= dims - 4; j += 4 )
                    {
                        float t0 = center[j] + w[i] * sample[j];
                        float t1 = center[j+1] + w[i] * sample[j+1];

                        center[j] = t0;
                        center[j+1] = t1;

                        t0 = center[j+2] + w[i] * sample[j+2];
                        t1 = center[j+3] + w[i] * sample[j+3];

                        center[j+2] = t0;
                        center[j+3] = t1;
                    }
                    for( ; j < dims; j++ )
                        center[j] += w[i] * sample[j];
                    counters[k]++;
                }

				//for( i = 0; i < N; i++ )
    //            {
    //                sample = data.ptr<float>(i);
    //                k = labels[i];
    //                float* center = centers.ptr<float>(k);
    //                for( j = 0; j < dims; j++ )
    //                {
				//		center[j] += w[i] * sample[j];
    //                }
    //                counters[k]++;
    //            }

                if( iter > 0 )
                    max_center_shift = 0;

                for( k = 0; k < K; k++ )
                {
                    float* center = centers.ptr<float>(k);
                    if( counters[k] != 0 )
                    {
                        float scale = 1.f/counters[k];
                        for( j = 0; j < dims; j++ )
                            center[j] *= scale;
                    }
                    else
                        generateRandomCenter(_box, center, rng);
                    
                    if( iter > 0 )
                    {
                        double dist = 0;
                        const float* old_center = old_centers.ptr<float>(k);
                        for( j = 0; j < dims; j++ )
                        {
                            double t = center[j] - old_center[j];
                            dist += t*t;
							//dist += fabs(t);
                        }
                        max_center_shift = std::max(max_center_shift, dist);
                    }
                }
            }

            // assign labels
            compactness = 0;
			//unsigned int *bins=new unsigned int[K];for( i = 0; i< K; i++) {bins[i]=0;}
            for( i = 0; i < N; i++ )
            {
                sample = data.ptr<float>(i);
                int k_best = 0;
                double min_dist = DBL_MAX;

                for( k = 0; k < K; k++ )
                {
                    const float* center = centers.ptr<float>(k);
                    double dist = mydistance(sample, center, dims);

                    if( min_dist > dist )
                    {
                        min_dist = dist;
                        k_best = k;
                    }
                }

                compactness += min_dist;
                labels[i] = k_best;
				//bins[k_best]++;
            }

			//for( k = 0; k< K; k++) {compactness+=abs(bins[k]-(float)N/K);}
			//delete []bins;
        }

        if( compactness < best_compactness )
        {
            best_compactness = compactness;
            if( _centers.needed() )
                centers.copyTo(_centers);
            _labels.copyTo(best_labels);
        }
    }

    return best_compactness;
}






int i,j;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	float *inData; 
	float *outData;
	float *weight;
	unsigned int *outLabel;
	unsigned int N=10000,C=10,D=128; 
	C = (unsigned int)mxGetScalar(prhs[1]);
	inData=(float *)mxGetPr(prhs[0]); 
	N=mxGetM(prhs[0]); 
	D=mxGetN(prhs[0]); 



	weight=(float *)mxGetPr(prhs[2]);


	cv::Mat FeatureMatrix(N,D,DataType<float>::type);
	for(i=0;i<N;i++){
		for(j=0;j<D;j++){ 
			FeatureMatrix.at<float>(i,j)=inData[j*N+i]; 
		}
	}
	delete [] prhs;
	FeatureMatrix.convertTo(FeatureMatrix,CV_32FC1,1,0);

	if( FeatureMatrix.empty()){
		cout<<"Input empty";
		return;
	}

	cv::Mat ClusterMatrix(C,D,DataType<float>::type);
	for( j = 0; j < C; j++ ){
	  for ( i = 0; i < D; i++ )
	  {
	   ClusterMatrix.at<float>(j,i)= 0.0f ;
	  }
	}


	cv::Mat label(N,1,DataType<unsigned int>::type);
	TermCriteria termcriteria(1,10,1e-6);

	mykmeans(FeatureMatrix,C,label,termcriteria,weight,1,KMEANS_PP_CENTERS ,ClusterMatrix);//KMEANS_PP_CENTERS  KMEANS_RANDOM_CENTERS



	
	plhs[0]=mxCreateNumericMatrix(C,D,mxSINGLE_CLASS,0); 
	outData=(float *)mxGetPr(plhs[0]); 
	for(i=0;i<C;i++){
		for(j=0;j<D;j++){ 
			outData[j*C+i]=ClusterMatrix.at<float>(i,j); 
		}
	}

	plhs[1]=mxCreateNumericMatrix(N,1,mxINT32_CLASS,0); 
	outLabel=(unsigned int *)mxGetPr(plhs[1]);
	for(i=0;i<N;i++){
		outLabel[i]=label.at<unsigned int>(i); 
	}

	FeatureMatrix.release();
	label.release();
	ClusterMatrix.release();
}