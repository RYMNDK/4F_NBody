#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <mpi.h>
#include <math.h>
#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"
#include <stdint.h>
#include <string.h>
#include <omp.h>

#define epsilon 0.000000000000022
#define G 6.67e-4
#define INDEX_POS 0
#define INDEX_V 1
#define INDEX_F 2
#define BORDER_REFLECT_X 0.1
#define BORDER_REFLECT_Y 0.1

#define _NUM_THREADS 4

#define __DEBUG 1

vec3* getItemsFrom2DArray(vec3* array, int index_X, int index_Y, int nRow, int nCol)
{
    
    if(index_X >= 0 && index_X < nRow){
        if(index_Y >= 0 && index_Y < nCol){
            return (array + index_X * nCol + index_Y);
        }
    }
    return nullptr;
}


void forceCalculation(vec3 * data, int offset, int numParticle, double * mass);
double relativeForce(vec3 * data, const double * mass, int particle1, int particle2, int numParticle ,int direction);
void fnPrintBuffer(void* buffer, int nAmount);
void fnPrintVEC(vec3* VEC, int nAmount);

void fnPaintPoint(char* img, long x, long y, int nWidth, int nHeight){
    
}

int world_rank;
int world_size;

double previousPositionX = -1.0;
double previousPositionY = -1.0;

int nImgIndex = 0;

double dStartTime = 0.0f;
double dEndTime = 0.0f;
double dTotalTime = 0.0f;

int main(int argc, char* argv[]){
    
    omp_set_num_threads(_NUM_THREADS);
    //this section initializes number of processors and current rank
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    

    if( argc != 10){
        printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
    }

    //variables
    int numParticlesLight = atoi(argv[1]);
    int numParticleMedium = atoi(argv[2]);
    int numParticlesHeavy = atoi(argv[3]);
    int numParticle = numParticlesLight + numParticleMedium + numParticlesHeavy;

    char* pCompleteFileName = (char*)calloc(sizeof(char), strlen(argv[9]) + 5 + 1 + 1 + 4);

    int numSteps = atoi(argv[4]);
    int subSteps = atoi(argv[5]);
    double timeSubStep = atof(argv[6]);

    int width = atoi(argv[7]);
    int height = atoi(argv[8]);

	double dWidth = width;
	double dHeight = height;

    auto * image = (unsigned char*)malloc(height*width*3*sizeof(unsigned char));

	int workChunk = numParticle / world_size;

	vec3 * receiveBuffer = (vec3 *)malloc(sizeof(vec3) * workChunk * 3);

	vec3 * sendBuffer = (vec3 *)malloc(sizeof(vec3) * workChunk * 3);

	
    
    
    //Initalize the image buffer.
    if (image != nullptr) {
        memset(image, 0x0, height*width * 3 * sizeof(unsigned char));
    }

    //auto * VEC = new vec3[numParticle*3];
    vec3* VEC = (vec3*)calloc(numParticle * 3, sizeof(vec3));
    double MASS[numParticle];
    //vec3 VEC[numParticle][3];
    //MASS[i] contains i's mass
    //VEC[i][0] contains i's position
    //VEC[i][1] contains i's velocity
    //VEC[i][2] contains the total force acting on i
    //initialize the 2d array of vec3
    if (world_rank == 0) {
        //this is a temporary vec3 to store a particle
        vec3 randVec;
        vec3 * tempVec;

        //randomize all the particles
        for (int i = 0; i < numParticle; i++){
            //for now each particle has random position and velocity
            //z is 0 for everything to make it easier

            //randomize the velocity of each particle
            if (i<numParticlesLight){
                MASS[i] = (massLightMax - massLightMin)*drand48() + massLightMin;
            
            }
            else if (i<numParticleMedium + numParticlesLight){
                MASS[i] = (massMediumMax - massMediumMin)*drand48() + massMediumMin;
                
            }
            else{
                MASS[i] = (massHeavyMax - massHeavyMin)*drand48() + massHeavyMin;
               
            }

            //randomize the position of each particle
            randVec = vec3((width - 1) * drand48(), (height - 1) * drand48(), 0);

            tempVec = getItemsFrom2DArray(VEC, i, INDEX_POS, numParticle, 3);
            if (tempVec != nullptr){
                *tempVec = randVec;

            } else {
                printf("initializing particle failed");
                if (VEC != nullptr) {
                    free(VEC);
                    VEC = nullptr;
                }
                if (image != nullptr) {
                    free(image);
                    image = nullptr;
                }
                MPI_Finalize();
                return -1;
            }
        }
    }

    //************** Finish Initialization of all the Particles ******************//


    //send out data to each processor
    MPI_Bcast(VEC, sizeof(vec3) * numParticle * 3, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(MASS, numParticle, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
	//printf("workChunk:%d, numParticle:%d, world_size:%d\r\n", workChunk, numParticle, world_size);

    for (int i = 0; i < numSteps * subSteps; i ++) {
		//printf("R:%d, 1\r\n", world_rank);
        //calculation
        
        if(world_rank == 0){
            dStartTime = MPI_Wtime();
        }
        #pragma omp parallel for
        for (int j = 0 + workChunk * world_rank; j < (world_rank + 1) * workChunk; j ++) {
			//printf("R:%d, 2\r\n", world_rank);
            //void forceCalculation(vec3 * data, int offset, int numParticle, double * mass)
            forceCalculation(VEC, j, numParticle, MASS);
            vec3 * currentForce = getItemsFrom2DArray(VEC, j, INDEX_F, numParticle, 3);
            vec3 * currentVelocity = getItemsFrom2DArray(VEC, j, INDEX_V, numParticle, 3);
            vec3 * currentPosition = getItemsFrom2DArray(VEC, j, INDEX_POS, numParticle, 3);
            double currentMass = MASS[j];
            double currentXa = currentForce->x / currentMass;
            double currentYa = currentForce->y / currentMass;
            //printf("currentForce-x: %f, currentForce-y: %f , currentXa:%f, currentYa:%f\r\n", currentForce->x, currentForce->y, currentXa, currentYa);
            
            //currentPosition->x += currentVelocity->x * timeSubStep;
            
            //currentPosition->y += currentVelocity->y * timeSubStep;
            
            
            
			currentVelocity->x += currentXa * timeSubStep;
			currentVelocity->y += currentYa * timeSubStep;

			double current_x = currentPosition->x;
			double current_y = currentPosition->y;

			//double Space_x = (currentVelocity->x * timeSubStep + currentXa * timeSubStep * timeSubStep *0.5);
			//double Space_y = (currentVelocity->y * timeSubStep + currentYa * timeSubStep * timeSubStep *0.5);
            double Space_x = currentVelocity->x * timeSubStep;
            double Space_y =  currentVelocity->y * timeSubStep;

            
            
			if (dWidth - current_x <= Space_x) {
				Space_x = fmod((Space_x - (dWidth - current_x)), dWidth);
				currentPosition->x = dWidth - Space_x;
				currentVelocity->x = -currentVelocity->x;
				//printf("1 x coor: %f, y coor: %f \r\n", currentPosition->x, currentPosition->y);
                currentVelocity->x *= BORDER_REFLECT_X;
                currentVelocity->y *= BORDER_REFLECT_Y;
			}
			else if (Space_x <= 0 - current_x) {
				Space_x = fmod((0 - current_x) - Space_x, dWidth);
				currentPosition->x = 0 + Space_x;
				currentVelocity->x = -currentVelocity->x;
				//printf("2 x coor: %f, y coor: %f \r\n", currentPosition->x, currentPosition->y);
                currentVelocity->x *= BORDER_REFLECT_X;
                currentVelocity->y *= BORDER_REFLECT_Y;
			}
			else
			{
				currentPosition->x += Space_x;
			}

			if (dHeight - current_y <= Space_y) {

				Space_y = fmod((Space_y - (dHeight - current_y)), dHeight);
				currentPosition->y = dHeight - Space_y;
				currentVelocity->y = -currentVelocity->y;
				//printf("3 x coor: %f, y coor: %f \r\n", currentPosition->x, currentPosition->y);
                currentVelocity->x *= BORDER_REFLECT_X;
                currentVelocity->y *= BORDER_REFLECT_Y;
			}
			else if (Space_y <= 0 - current_y) {
				Space_y = fmod((0 - current_y) - Space_y, dHeight);
				currentPosition->y = 0 + Space_y;
				currentVelocity->y = -currentVelocity->y;
				//printf("4 x coor: %f, y coor: %f \r\n", currentPosition->x, currentPosition->y);
                currentVelocity->x *= BORDER_REFLECT_X;
                currentVelocity->y *= BORDER_REFLECT_Y;
			}
			else
			{
				currentPosition->y += Space_y;
			}
             
            
			

			/*printf("index:%d, mass:%f, theta_x:%f, theta_y: %f, x:%f, y:%f\r\n", 
				j,  
				MASS[j],
				Space_x,
				Space_y,
				currentPosition->x,
				currentPosition->y);**/
            //printf("x coor: %d, y coor: %d \r\n", currentPosition->x, currentPosition->y);

			




            /*

            if (currentPosition->x >= (double)(width - 1)) {
                currentPosition->x = fmod(currentPosition->x, width - 1);
                currentVelocity->x *= 0.8;
                
                
				//Reflection

				//printf("1.1 x:%f, y:%f\r\n", currentPosition->x, currentPosition->y);
				
            }else if (currentPosition->x <= 1) {
                currentPosition->x = ((width - 1)  - fmod(abs(currentPosition->x), (width - 1)));
                //printf("2 x:%f, y:%f\r\n", currentPosition->x, currentPosition->y);
                currentVelocity->x *= 0.8;
            }

            if (currentPosition->y >= (double)(height - 1)) {
                currentPosition->y = fmod(currentPosition->y, height - 1);
                //if(nImgIndex == 16){
               //     printf("3 x:%f, y:%f, height:%f\r\n", currentPosition->x, currentPosition->y, (double)height);                }
                currentVelocity->y *= 0.8;
                
            }else if (currentPosition->y <= 1) {
                currentPosition->y = ((height - 1) - fmod(abs(currentPosition->y), (height - 1)));
                //printf("4 x:%f, y:%f\r\n", currentPosition->x, currentPosition->y);
                currentVelocity->y *= 0.8;
            }
             */
            
            //reset force

			
			currentForce->x = 0.0;
			currentForce->y = 0.0;
        }

        /*vec3 * firstParticlePosition = getItemsFrom2DArray(VEC, 0, INDEX_POS, numParticle, 3);

        if (firstParticlePosition->x == previousPositionX && firstParticlePosition->y == previousPositionY) {
            printf("the particle's position is not changed at all");
            return -1;
        } else {
            previousPositionX = firstParticlePosition->x;
            previousPositionY = firstParticlePosition->y;
        }**/

        
		memset(receiveBuffer, 0x0, sizeof(vec3) * workChunk * 3);
        
		memset(sendBuffer, 0x0, sizeof(vec3) * workChunk * 3);
        int nChrunkIndex = 0;
        
        
        MPI_Request reqs[4];
        MPI_Status stats[4];
        //for world_size-1 times communication

        //get the address of the starting element of a chunk of a data and use offset to fetch all the chunk
        memcpy(
                sendBuffer,
               getItemsFrom2DArray(VEC, workChunk * world_rank, INDEX_POS, numParticle, 3),
               sizeof(vec3)*workChunk*3);

        nChrunkIndex = world_rank;
        int nSendRevSize = (sizeof(vec3)/sizeof(double)) * workChunk * 3;
        
        

        for (int k = 0; k < world_size -1 ; k++) {
			
            //each processor has to send out data it holds in each iteration
            MPI_Isend(sendBuffer,
                          nSendRevSize,
                          MPI_DOUBLE,
                          (world_rank+1)%world_size,
                          world_rank,
                          MPI_COMM_WORLD,
                          &reqs[0]);
           
            MPI_Isend(&nChrunkIndex,
                      1,
                      MPI_INT,
                      (world_rank+1)%world_size,
                      world_rank + 10000,
                      MPI_COMM_WORLD,
                      &reqs[1]);
            


            MPI_Irecv(&nChrunkIndex,
                      1,
                      MPI_INT,
                      (world_rank + world_size - 1) % world_size,
                      (world_rank + world_size - 1) % world_size + 10000,
                      MPI_COMM_WORLD,
                      &reqs[2]);
 
            
            MPI_Irecv(receiveBuffer,
                          nSendRevSize,
                          MPI_DOUBLE,
                          (world_rank + world_size - 1) % world_size,
                          (world_rank + world_size - 1) % world_size,
                          MPI_COMM_WORLD,
                          &reqs[3]);
            

            //printf("Rank:%d, receivefrom:%d, sendto:%d\r\n", world_rank, (world_rank + world_size - 1) % world_size, (world_rank+1)%world_size);

            MPI_Waitall(4, reqs, stats);
            
            //Update local buffer.

            
            memcpy(
				VEC + nChrunkIndex  * workChunk * 3,
				receiveBuffer,
				sizeof(vec3) * workChunk * 3
			);

            //Prepare next round send buffer.
            memcpy(sendBuffer, receiveBuffer, sizeof(vec3) * workChunk * 3);
            
	/*if(nImgIndex == 15){
		printf("-----------\r\n");
		fnPrintVEC(receiveBuffer, workChunk);
		printf("***********\r\n");
		fnPrintVEC(VEC, numParticle);
	} **/     
        


        }
        
        if(world_rank == 0){
            dEndTime = MPI_Wtime();
            dTotalTime += dEndTime - dStartTime;
        }

        
        
        
        
        //record state
        if (world_rank == 0 && i % subSteps == 0) {
            nImgIndex = i / subSteps;
			//printf("R:%d, 4\r\n", world_rank);
            if (image != nullptr) {
                memset(image, 0xFF, height * width * 3 * sizeof(unsigned char));
            }
            
            
            
            //fnPrintVEC(VEC, numParticle);
            #pragma omp parallel for
            for (int i = 0; i < numParticle; i++) {
                vec3* position_ = getItemsFrom2DArray(VEC, i, INDEX_POS, numParticle, 3);
                int x_ = position_->x;
                int y_ = position_->y;

                //printf("particle_index:%d, x:%f, y:%f\r\n", i, position_->x, position_->y);

                //int nPosition = (VEC[i * 3 + 0].y * width + VEC[i * 3 + 0].x) * 3;
                //int nPosition = (y_ * width + x_) * 3;
                
				/*
                if(nImgIndex == 2){
                    printf("nPosition: %d, i:%d, x:%d, y:%d\r\n", nPosition, i, x_, y_);
                }*/
                

                if (MASS[i] >= massLightMin && MASS[i] <= massLightMax) {
                    //printf("Light\r\n");
                    //auto * image = (unsigned char*)malloc(height*width*3*sizeof(unsigned char));
                    
                    if (y_  > 0.0 && y_ < (height - 1) && x_> 0.0 && x_ < (width - 1)) {
                        //left to the current
                        image[(y_ * width + x_ - 1) * 3 + 0] = 0xFF;
                        image[(y_ * width + x_ - 1 ) * 3 + 1] = 0x0;
                        image[(y_ * width + x_ - 1) * 3 + 2] = 0x0;
                        
                        //right to the current
                        image[(y_ * width + x_ + 1) * 3 + 0] = 0xFF;
                        image[(y_ * width + x_ + 1 ) * 3 + 1] = 0x0;
                        image[(y_ * width + x_ + 1) * 3 + 2] = 0x0;
                        
                        //above the current
                        image[((y_ -1) * width + x_) * 3 + 0] = 0xFF;
                        image[((y_ -1)* width + x_) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_) * 3 + 2] = 0x0;
                        
                        //below the current
                        image[((y_ +1) * width + x_) * 3 + 0] = 0xFF;
                        image[((y_ +1)* width + x_) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_) * 3 + 2] = 0x0;
                        
                        //left to the above
                        image[((y_ -1) * width + x_ -1) * 3 + 0] = 0xFF;
                        image[((y_ -1)* width + x_ -1) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_-1) * 3 + 2] = 0x0;
                        
                        //right to the above
                        image[((y_ -1) * width + x_+1) * 3 + 0] = 0xFF;
                        image[((y_ -1)* width + x_+1) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_+1) * 3 + 2] = 0x0;
                        
                        //left to the below
                        image[((y_ +1) * width + x_ -1) * 3 + 0] = 0xFF;
                        image[((y_ +1)* width + x_ -1) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_-1) * 3 + 2] = 0x0;
                        
                        //right to the below
                        image[((y_ +1) * width + x_+1) * 3 + 0] = 0xFF;
                        image[((y_ +1)* width + x_+1) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_+1) * 3 + 2] = 0x0;
                    }
                    
                    //current particle
                    image[(y_ * width + x_) * 3 + 0] = 0xFF;
                    image[(y_ * width + x_) * 3 + 1] = 0x0;
                    image[(y_ * width + x_) * 3 + 2] = 0x0;
                
                    

                }
                else if (MASS[i] >= massMediumMin && MASS[i] <= massMediumMax) {
                    if (y_ - 1  >= 0.0 && y_ < height - 1 && x_-1 >= 0.0 && x_ < width - 1){
                        //left to the current
                        image[(y_ * width + x_ - 1) * 3 + 0] = 0x0;
                        image[(y_ * width + x_ - 1 ) * 3 + 1] = 0xFF;
                        image[(y_ * width + x_ - 1) * 3 + 2] = 0x0;
                        
                        //right to the current
                        image[(y_ * width + x_ + 1) * 3 + 0] = 0x0;
                        image[(y_ * width + x_ + 1 ) * 3 + 1] = 0xFF;
                        image[(y_ * width + x_ + 1) * 3 + 2] = 0x0;
                        
                        //above the current
                        image[((y_ -1) * width + x_) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_) * 3 + 1] = 0xFF;
                        image[((y_-1) * width + x_) * 3 + 2] = 0x0;
                        
                        //below the current
                        image[((y_ +1) * width + x_) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_) * 3 + 1] = 0xFF;
                        image[((y_+1) * width + x_) * 3 + 2] = 0x0;
                        
                        //left to the above
                        image[((y_ -1) * width + x_ -1) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_ -1) * 3 + 1] = 0xFF;
                        image[((y_-1) * width + x_-1) * 3 + 2] = 0x0;
                        
                        //right to the above
                        image[((y_ -1) * width + x_+1) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_+1) * 3 + 1] = 0xFF;
                        image[((y_-1) * width + x_+1) * 3 + 2] = 0x0;
                        
                        //left to the below
                        image[((y_ +1) * width + x_ -1) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_ -1) * 3 + 1] = 0xFF;
                        image[((y_+1) * width + x_-1) * 3 + 2] = 0x0;
                        
                        //right to the below
                        image[((y_ +1) * width + x_+1) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_+1) * 3 + 1] = 0xFF;
                        image[((y_+1) * width + x_+1) * 3 + 2] = 0x0;
                    }
                    //printf("Medium\r\n");
                    image[(y_ * width + x_) * 3 + 1] = 0xFF;
                    image[(y_ * width + x_) * 3 + 0] = 0x0;
                    image[(y_ * width + x_) * 3 + 2] = 0x0;
                    //printf("nPosition: %ld, i:%d, x:%ld, y:%ld\r\n", nPosition, i, x_, y_);
                }
                else if (MASS[i] >= massHeavyMin && MASS[i] <= massHeavyMax) {
                    if (y_ - 1  >= 0.0 && y_ < height - 1 && x_-1 >= 0.0 && x_ < width - 1){
                        //left to the current
                        image[(y_ * width + x_ - 1) * 3 + 0] = 0x0;
                        image[(y_ * width + x_ - 1 ) * 3 + 1] = 0x0;
                        image[(y_ * width + x_ - 1) * 3 + 2] = 0xFF;
                        
                        //right to the current
                        image[(y_ * width + x_ + 1) * 3 + 0] = 0x0;
                        image[(y_ * width + x_ + 1 ) * 3 + 1] = 0x0;
                        image[(y_ * width + x_ + 1) * 3 + 2] = 0xFF;
                        
                        //above the current
                        image[((y_ -1) * width + x_) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_) * 3 + 2] = 0xFF;
                        
                        //below the current
                        image[((y_ +1) * width + x_) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_) * 3 + 2] = 0xFF;
                        
                        //left to the above
                        image[((y_ -1) * width + x_ -1) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_ -1) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_-1) * 3 + 2] = 0xFF;
                        
                        //right to the above
                        image[((y_ -1) * width + x_+1) * 3 + 0] = 0x0;
                        image[((y_ -1)* width + x_+1) * 3 + 1] = 0x0;
                        image[((y_-1) * width + x_+1) * 3 + 2] = 0xFF;
                        
                        //left to the below
                        image[((y_ +1) * width + x_ -1) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_ -1) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_-1) * 3 + 2] = 0xFF;
                        
                        //right to the below
                        image[((y_ +1) * width + x_+1) * 3 + 0] = 0x0;
                        image[((y_ +1)* width + x_+1) * 3 + 1] = 0x0;
                        image[((y_+1) * width + x_+1) * 3 + 2] = 0xFF;
                    }
                    //printf("Heavy\r\n");
                    image[(y_ * width + x_) * 3 + 2] = 0xFF;
                    image[(y_ * width + x_) * 3 + 0] = 0x0;
                    image[(y_ * width + x_) * 3 + 1] = 0x0;
                }
            }
            
            pCompleteFileName = (char*)memset(pCompleteFileName, 0x0, sizeof(char) * (strlen(argv[9]) + 5 + 1 + 1 + 4) );
            sprintf(pCompleteFileName, "%s_%05d.bmp", argv[9], nImgIndex);
            saveBMP(pCompleteFileName, image, width, height);

        }

		

    }
    
    if(world_rank == 0){
        printf("numParticle:%d, world_size:%d, total time:%f\r\n", numParticle, world_size, dTotalTime);
    }
    


    if(pCompleteFileName != nullptr){
        free(pCompleteFileName);
        pCompleteFileName = nullptr;
    }
    
    free(image);
    free(VEC);
	free(receiveBuffer);
	free(sendBuffer);
    MPI_Finalize();
    return 0;
}

void forceCalculation(vec3 * data, int offset, int numParticle, double * mass) {
    vec3 * currentForce = getItemsFrom2DArray(data, offset, INDEX_F, numParticle, 3);
#pragma omp parallel for
    for (int i = offset + 1; i < numParticle; i++) {
        double x_relativeForce = relativeForce(data, mass, offset, i, numParticle, 0);
        double y_relativeForce = relativeForce(data, mass, offset, i, numParticle, 1);
        

        
        
        /*
        if(x_relativeForce <= 0)
        {
            printf("fx:%f\r\n", x_relativeForce);
        }
        
        if(y_relativeForce <= 0){
            printf("fy:%f\r\n", y_relativeForce);
        }*/

		/*if (nImgIndex == 19) {
			printf("x_r:%f, y_r:%f\r\n", x_relativeForce, y_relativeForce);
		}*/
        
        
        
        
        currentForce->x += x_relativeForce;
        currentForce->y += y_relativeForce;
        

        vec3 * otherForce = getItemsFrom2DArray(data, i, INDEX_F, numParticle, 3);
        otherForce->x -= x_relativeForce;
        otherForce->y -= y_relativeForce;
    }
}

inline double relativeForce(vec3 * data, const double * mass, int particle1, int particle2, int numParticle ,int direction) {
    double Mass1 = mass[particle1];
    double Mass2 = mass[particle2];

    vec3 * particle1direction = getItemsFrom2DArray(data, particle1, INDEX_POS, numParticle, 3);
    vec3 * particle2direction = getItemsFrom2DArray(data, particle2, INDEX_POS, numParticle, 3);

    double particle1X = particle1direction->x;
    double particle1Y = particle1direction->y;
    double particle2X = particle2direction->x;
    double particle2Y = particle2direction->y;

    double distance = sqrt((particle1X - particle2X)*(particle1X - particle2X)
                       +
                       (particle1Y - particle2Y)*(particle1Y - particle2Y));
	

    if (direction == 0) {
        double cos = (particle1X - particle2X) / distance;
        //printf("cos:%f\r\n",cos);
        return (G * Mass1 * Mass2 / (distance * distance  + epsilon * epsilon  ))  * cos;
    } else {
        double sin = (particle1Y - particle2Y) / distance;
        return (G * Mass1 * Mass2 / (distance * distance * distance + epsilon  * epsilon  ))  * sin;
    }
}


void fnPrintBuffer(void* buffer, int nAmount){
    if(buffer != NULL){
        printf("Rank:%d\r\n", *(int*)buffer);
        buffer = (char*)buffer + sizeof(int);
        for(int i = 0; i < nAmount; i++){
            printf("index:%d, x:%f, y:%f\r\n", i, (((vec3*)buffer) + i)->x, (((vec3*)buffer) + i)->y);
        }
    }
}

void fnPrintVEC(vec3* VEC, int nAmount){
    if(VEC != NULL){
        for(int i = 0; i < nAmount * 3; i+= 3){
            printf("index:%d, x: %f, y: %f, vx:%f, vy:%f, fx:%f, fy:%f\r\n",
                   i/3,
                   VEC[i + 0].x,
                   VEC[i + 0].y,
                   VEC[i + 1].x,
                   VEC[i + 1].y,
                   VEC[i + 2].x,
                   VEC[i + 2].y
            );
        }
    }
}
