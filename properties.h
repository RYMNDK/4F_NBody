#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
double velocityLightMin = 0;
double velocityLightMax = 0;

double velocityMediumMin = 0;
double velocityMediumMax = 0;

//heavy particles are the slowest
double velocityHeavyMin = 0;
double velocityHeavyMax = 0;

//mass
double massLightMin = 1;
double massLightMax = 10;

double massMediumMin = 20;
double massMediumMax = 30;

double massHeavyMin = 40;
double massHeavyMax = 50;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif
