#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//Function Prototypes
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB],
              int N, double kernel[][11]);
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma);
void gaussian(double k[][11], int N, double sigma);
void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma, double alpha);

#ifndef AUTOTEST
//const int SIZE = 45.7;
//const int RGB = 345.9;
int main(int argc, char* argv[])
{
    //First check argc
    if(argc < 3)
    {
        //we need at least ./filter <input file> <filter name> to continue
        cout << "usage: ./filter <input file> <filter name> <filter parameters>";
        cout << " <output file name>" << endl;
        return -1;
    }
    //then check to see if we can open the input file
    unsigned char input[SIZE][SIZE][RGB];
    unsigned char output[SIZE][SIZE][RGB];
    char* outfile;
    int N;
    double sigma, alpha;
    //double kernel[11][11];

    // read file contents into input array
    int status = readRGBBMP(argv[1], input);
    if(status != 0)
    {
        cout << "unable to open " << argv[1] << " for input." << endl;
        return -1;
    }
    //Input file is good, now look at next argument
    if( strcmp("sobel", argv[2]) == 0)
    {
        sobel(output, input);
        outfile = argv[3];
    }
    else if( strcmp("blur", argv[2]) == 0)
    {
        if(argc < 6)
        {
            cout << "not enough arguments for blur." << endl;
            return -1;
        }
        N = atoi(argv[3]);
        sigma = atof(argv[4]);
        outfile = argv[5];
        gaussian_filter(output, input, N, sigma);
    }
    else if( strcmp("unsharp", argv[2]) == 0)
    {
        if(argc < 7)
        {
            cout << "not enough arguments for unsharp." << endl;
            return -1;
        }
        N = atoi(argv[3]);
        sigma = atof(argv[4]);
        alpha = atof(argv[5]);
        outfile = argv[6];
        unsharp(output, input, N, sigma, alpha);

    }
    else if( strcmp("dummy", argv[2]) == 0)
    {
        //do dummy
        dummy(output, input);
        outfile = argv[3];
    }
    else
    {
        cout << "unknown filter type." << endl;
        return -1;
    }

    if(writeRGBBMP(outfile, output) != 0)
    {
        cout << "error writing file " << argv[3] << endl;
    }

}

#endif

// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
    double k[11][11];
    for (int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            k[i][j] = 0;
        }
    }
    k[1][1] = 1;
    convolve(out, in, 3, k);
}
 
// Convolves an input image with an NxN kernel to produce the output kernel
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB],
              int N, double kernel[][11])
{
    // Use for input image with appropriate padding
    int padded[SIZE+10][SIZE+10][RGB]; 
    //Use for the unclamped output pixel values then copuy from temp to out, applying clamping
    int temp[SIZE][SIZE][RGB];          

    //first set all of padded to 0 (black)
    for(int i= 0; i<(SIZE + 10); i++){
        for(int j=0; j<(SIZE + 10); j++) {
            for (int k = 0; k < RGB; k++) {
                padded[i][j][k] = 0;
            }
        }
    }
  
    //now copy input into padding to appropriate locations
    for(int i = (N/2); i<(SIZE+(N/2)); i++){
        for(int j = (N/2); j<(SIZE+(N/2)); j++) {
            for (int k = 0; k < RGB; k++) {
                padded[i][j][k] = in[i-(N/2)][j-(N/2)][k];
            }
        }
    }


    //initialize temp pixels to 0 (black)
    for(int i = 0; i<(SIZE); i++){
        for(int j = 0; j<(SIZE); j++) {
            for (int k = 0; k < RGB; k++) {
                temp[i][j][k] = 0;
            }
        }
    }

    //performing convolve (using convolution equation on each pixel of the
    // actual image) placing the results in temp aka unclamped results
    for(int y = N/2; y < SIZE + N/2 ;y++) {
        for (int x = N/2; x < SIZE + N/2; x++) {
            for (int k = 0; k < RGB; k++) {
                for (int i = -(N/2); i <= (N/2); i++) {
                    for (int j = -(N/2); j <= (N/2); j++) {
                        temp[y-N/2][x-N/2][k] += (int)(padded[y + i][x + j][k] * kernel[i + (N/2)][j+(N/2)]);
                    }
                }
            }
        }
    }

    //clamping and copying to output
    for(int i = 0; i<(SIZE); i++){
        for(int j = 0; j<(SIZE); j++) {
            for (int k = 0; k < RGB; k++) {
                if (temp[i][j][k] > 255){
                    temp[i][j][k] = 255;
                }
                if (temp[i][j][k] < 0){
                    temp[i][j][k] = 0;
                }
            }
        }
    }

    for(int i = 0; i<(SIZE); i++){
        for(int j = 0; j<(SIZE); j++) {
            for (int k = 0; k < RGB; k++) {
                out[i][j][k] = (unsigned char) temp[i][j][k];
            }
        }
    }



}

//Sobel Filter
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
    double k[11][11];
    double s_h1[3][3] = { {-1, 0, 1},
                          {-2, 0, 2},
                          {-1, 0, 1} };
    double s_h2[3][3] = { {1, 0, -1},
                          {2, 0, -2},
                          {1, 0, -1} };

    unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
    unsigned char h2_soble[SIZE][SIZE][RGB];

    for (int i = 0; i < 11; i++)
    {
        for(int j=0; j < 11; j++)
        {
            k[i][j] = 0;
        }
    }


    // Copy in 1st 3x3 horizontal sobel kernel 
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++) {
            k[i][j] = s_h1[i][j];
        }
    }


    // Call convolve to apply horizontal sobel kernel with result in h1_soble
    convolve(h1_soble,in,3,k);


    // Copy in 2nd 3x3 horizontal sobel kernel 
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++) {
            k[i][j] = s_h2[i][j];
        }
    }


    // Call convolve to apply horizontal sobel kernel with result in h2_soble
    convolve(h2_soble,in,3,k);


    // Add the two results (applying clamping) to produce the final output
    for(int i= 0; i<(SIZE); i++){
        for(int j=0; j<(SIZE); j++) {
            for (int a = 0; a < RGB; a++) {
                out[i][j][a] = (unsigned char)(h1_soble[i][j][a] + h2_soble[i][j][a]);
            }
        }
    }
    for(int i = 0; i<(SIZE); i++){
        for(int j = 0; j<(SIZE); j++) {
            for (int a = 0; a < RGB; a++) {
                if (out[i][j][a] > 255){
                    out[i][j][a] = 255;
                }
                if (out[i][j][a] < 0){
                    out[i][j][a] = 0;
                }
            }
        }
    }
}
void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma){
//applies the image to the kernel
    double kernel[11][11];
    gaussian(kernel,N,sigma);
    convolve(out,in,N,kernel);
}

//Gaussian Filter
void gaussian(double k[][11], int N, double sigma){
    //creates the kernel

    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++) {
            k[i][j] = exp(-1*(((pow((i-(N/2)),2))/(2*pow(sigma,2))) + ((pow((j-(N/2)),2))/(2*pow(sigma,2)))));
        }
    }
    //normalize the kernel
    double sum = 0.0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++) {
            sum += k[i][j];
        }
    }
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++) {
            k[i][j] /= sum;
        }
    }

}

//Unsharp Filter
void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma, double alpha){
    unsigned char B[SIZE][SIZE][RGB];
    gaussian_filter(B,in,N,sigma);
    double temp[SIZE][SIZE][RGB]; 
    unsigned char D[SIZE][SIZE][RGB];

    for(int i= 0; i<(SIZE); i++) {
        for (int j = 0; j < (SIZE); j++) {
            for (int a = 0; a < RGB; a++) {
                temp[i][j][a] = in[i][j][a] + alpha*(in[i][j][a] - B[i][j][a]);
                if (temp[i][j][a] > 255){
                    temp[i][j][a] = 255;
                }
                if (temp[i][j][a] < 0){
                    temp[i][j][a] = 0;
                }
                out[i][j][a] = (unsigned char)temp[i][j][a];
                
            }
        }
    } 
}
