//************************************************************
// Demo OpenCL application to compute a simple vector addition
// computation between 2 arrays on the GPU
// ************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
#include <time.h>
// OpenCL source code
const char* OpenCLSource = "\n" \
"__kernel void VectorAdd(__global int* c, __global int* a,__global int* b,"\
" const unsigned int nmax)"\
"{\n"\
"      unsigned int n = get_global_id(0);\n"\
"      if( n < nmax ){ c[n] = a[n] + b[n]; } \n"\
"}\n";


// Some interesting data for the vectors
int InitialData1[20] = {37,50,54,50,56,0,43,43,74,71,32,36,16,43,56,100,50,25,15,17
};
int InitialData2[20] = {35,51,54,58,55,32,36,69,27,39,35,40,16,44,55,14,58,75,18,15
};

// Number of elements in the vectors to be added
#define SIZE 2048*2048*5

void die(char * msg){
    puts(msg);
    exit(0);
}

#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define CLDO(f) (f == CL_SUCCESS ) ? (void) 0 \
                    : die("OpenCL Error in "__FILE__":"_STR(__LINE__) " " #f  ) 

// Main function
// ************************************************************
int main(int argc, char **argv)
{
  // Two integer source vectors in Host memory
  int i, c;
  cl_platform_id *cpPlatform;
  cl_device_id cdDevice;
  char cBuffer[1024];
  cl_context GPUContext;
  cl_command_queue cqCommandQueue;
  cl_mem GPUVector1 , GPUVector2, GPUOutputVector ;
  cl_program OpenCLProgram;
  cl_kernel OpenCLVectorAdd;
  cl_int err, count;
  cl_uint num_platforms = 0, num_devices = 0;
  size_t WorkSize[1] = {SIZE};
  int *HostVector1, *HostVector2;
  int *HostOutputVector, *testVector;

  clock_t start, end;

  start = clock();
  count = SIZE;
  // Initialize with some interesting repeating data
  HostVector1 = malloc( SIZE * sizeof(int) );
  HostVector2 = malloc( SIZE * sizeof(int) );
  HostOutputVector = malloc( SIZE * sizeof(int) );
  testVector  = malloc( SIZE * sizeof(int) );
  
  c = 0;
  while( c < SIZE )
  {
         HostVector1[c] = InitialData1[c%20];
         HostVector2[c] = InitialData2[c%20];
         testVector[c] = HostVector1[c]+HostVector2[c];
         c++;
  }

  end = clock();
  printf("Init and CPU %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();

  // Find out how many opencl plaforms are available
  CLDO( clGetPlatformIDs( 0, NULL, &num_platforms));
  printf("Found %d opencl platforms\n", num_platforms);

  cpPlatform = malloc( num_platforms * sizeof(cl_platform_id) );
  if ( cpPlatform == NULL ){
      printf("Malloc failed\n");
      exit(0);
  }
  
  //Get an OpenCL platform
  CLDO( clGetPlatformIDs(2, cpPlatform, NULL) );

  // List OpenCL devices
//  for( i=num_platforms-1 ; i >=0 ; i--) { 
  for( i=0 ; i < num_platforms ; i++) { 
     CLDO(clGetDeviceIDs(cpPlatform[i], CL_DEVICE_TYPE_ALL, 1, &cdDevice, &num_devices));
     printf("Platform %d has %d devices\n",i,num_devices);
     CLDO(clGetDeviceInfo(cdDevice, CL_DEVICE_NAME, sizeof(cBuffer), &cBuffer, NULL));
     printf("CL_DEVICE_NAME:       %s\n", cBuffer);
     CLDO(clGetDeviceInfo(cdDevice, CL_DRIVER_VERSION, sizeof(cBuffer), &cBuffer, NULL));
     printf("CL_DRIVER_VERSION: %s\n\n", cBuffer);
  }


  // Create a context to run OpenCL on our CUDA-enabled NVIDIA GPU
  GPUContext = clCreateContext(0, 1, &cdDevice, NULL, NULL, &err);
  CLDO(err);
  
  // Create a command-queue on the GPU device
  cqCommandQueue = clCreateCommandQueue(GPUContext, cdDevice, 0, &err);
  CLDO(err);

  end = clock();
  printf("Init GPU %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();


  // Allocate GPU memory for source vectors AND initialize from CPU memory
  GPUVector1 = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY |
                                CL_MEM_COPY_HOST_PTR, sizeof(int) * SIZE, HostVector1, &err); 
  CLDO(err);
  GPUVector2 = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY |
                                 CL_MEM_COPY_HOST_PTR, sizeof(int) * SIZE, HostVector2, &err); 
  CLDO(err);

  // Allocate output memory on GPU
  GPUOutputVector = clCreateBuffer(GPUContext, CL_MEM_WRITE_ONLY,
                                       sizeof(int) * SIZE, NULL, &err); 
  CLDO(err);
  end=clock();
  printf("Allocated buffers on device\n",
          ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();
  // Create OpenCL program with source code
  //
  printf("Going to try to compile\nProgram source is:\n%s", OpenCLSource);

  OpenCLProgram = clCreateProgramWithSource(GPUContext, 1,
                                  &OpenCLSource, NULL, &err); 
  CLDO(err);

  // Build the program (OpenCL JIT compilation)
  err = clBuildProgram(OpenCLProgram, 0, NULL, NULL, NULL, NULL);
  if( err ) {
      printf("Compile failure.\n");
      exit(0);
  }

  // Create a handle to the compiled OpenCL function (Kernel)
  OpenCLVectorAdd = clCreateKernel(OpenCLProgram, "VectorAdd", &err); CLDO(err);

  end = clock();
  printf("Alloc & compile %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();


  // In the next step we associate the GPU memory with the Kernel arguments
  CLDO(clSetKernelArg(OpenCLVectorAdd, 0, sizeof(cl_mem), (void*)&GPUOutputVector) );
  CLDO(clSetKernelArg(OpenCLVectorAdd, 1, sizeof(cl_mem), (void*)&GPUVector1));
  CLDO(clSetKernelArg(OpenCLVectorAdd, 2, sizeof(cl_mem), (void*)&GPUVector2));
  CLDO(clSetKernelArg(OpenCLVectorAdd, 3, sizeof(cl_uint), (void*)&count ));

  end = clock();
  printf("Copy data %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();

  // Launch the Kernel on the GPU
  // one dimensional Range
  CLDO(clEnqueueNDRangeKernel(cqCommandQueue, OpenCLVectorAdd, 1, NULL,
                              WorkSize, NULL, 0, NULL, NULL));

  end = clock();
  printf("Run kernel %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();


  // Copy the output in GPU memory back to CPU memory
  CLDO( clEnqueueReadBuffer(cqCommandQueue, GPUOutputVector, CL_TRUE, 0,
            SIZE * sizeof(int), HostOutputVector, 0, NULL, NULL));

  end = clock();
  printf("Read GPU %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();

  // Cleanup
  CLDO(clReleaseKernel(OpenCLVectorAdd));
  CLDO(clReleaseProgram(OpenCLProgram));
  CLDO(clReleaseCommandQueue(cqCommandQueue));
  CLDO(clReleaseContext(GPUContext));
  CLDO(clReleaseMemObject(GPUVector1));
  CLDO(clReleaseMemObject(GPUVector2));
  CLDO(clReleaseMemObject(GPUOutputVector));
  free(cpPlatform);

  end = clock();
  printf("Free GPU %f\n", ( (double) end - (double) start ) / (double) CLOCKS_PER_SEC);
  start = clock();

  // Print out the results
  for(c = 0; c <SIZE ; c++)
    if( HostOutputVector[ c ] != testVector[c] )
        printf("Error %d\n",c);
  


  printf("\n\nThe End\n\n");
  return 0;
}

