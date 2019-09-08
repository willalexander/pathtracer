#include <random>
#include <chrono>

#include "randomnumbers.h"

float randNumbers[2000];
int randCounter;

static int rCount = 0;

std::default_random_engine generator(rCount);
//rCount++;

void generateRandomNumbers()
{
	static int rCount = 0;
	
	std::default_random_engine generator(rCount);
	rCount++;
	std::uniform_real_distribution<float> distribution(0.0, 1.0);
	
	for(int i = 0; i < 2000; i++)
	{
		randNumbers[i] = distribution(generator);
	}
	
	randCounter = 0;
}

float randFloat()
{
	/*randCounter++;
	if(randCounter > 1000) randCounter = 0;
	return randNumbers[randCounter];*/

	
	
	std::uniform_real_distribution<float> distribution(0.0, 1.0);

	return distribution(generator);
}

bool randBool()
{
	if (randFloat() < 0.5) return true;
	return false;
}