#define cufftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "cufft.h"
#include "stdio.h"
    static const char *_hipGetErrorEnum(cufftResult error)
    {
    switch (error)
    {
    case CUFFT_SUCCESS:
    return "CUFFT_SUCCESS";

    case CUFFT_INVALID_PLAN:
    return "CUFFT_INVALID_PLAN";

    case CUFFT_ALLOC_FAILED:
    return "CUFFT_ALLOC_FAILED";

    case CUFFT_INVALID_TYPE:
    return "CUFFT_INVALID_TYPE";

    case CUFFT_INVALID_VALUE:
    return "CUFFT_INVALID_VALUE";

    case CUFFT_INTERNAL_ERROR:
    return "CUFFT_INTERNAL_ERROR";

    case CUFFT_EXEC_FAILED:
    return "CUFFT_EXEC_FAILED";

    case CUFFT_SETUP_FAILED:
    return "CUFFT_SETUP_FAILED";

    case CUFFT_INVALID_SIZE:
    return "CUFFT_INVALID_SIZE";

    case CUFFT_UNALIGNED_DATA:
    return "CUFFT_UNALIGNED_DATA";

    case CUFFT_INCOMPLETE_PARAMETER_LIST:
    return "CUFFT_INCOMPLETE_PARAMETER_LIST";

    case CUFFT_INVALID_DEVICE:
    return "CUFFT_INVALID_DEVICE";

    case CUFFT_PARSE_ERROR:
    return "CUFFT_PARSE_ERROR";

    case CUFFT_NO_WORKSPACE:
    return "CUFFT_NO_WORKSPACE";

    case CUFFT_NOT_IMPLEMENTED:
    return "CUFFT_NOT_IMPLEMENTED";

    case CUFFT_NOT_SUPPORTED:
    return "CUFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
    }

    inline void __hipfftSafeCall(cufftResult err, const char *file, const int line)
    {
    if( CUFFT_SUCCESS != err) {
    fprintf(stderr, "CUFFT error at 1\n");
    fprintf(stderr, "CUFFT error in file '%s'\n",__FILE__);
    fprintf(stderr, "CUFFT error at 2\n");
    /*fprintf(stderr, "CUFFT error line '%s'\n",__LINE__);*/
    fprintf(stderr, "CUFFT error at 3\n");
    /*fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
    _hipGetErrorEnum(err)); \*/
    fprintf(stderr, "CUFFT error %d: %s\nterminating!\n",err,_hipGetErrorEnum(err)); \
    cudaDeviceReset(); return; \
    }
    }

extern "C"
void
destroy_plan_ffth_(cufftHandle *PLANp)
{
cufftHandle plan = *PLANp;

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

cufftSafeCall(cufftDestroy(plan));

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}


}

