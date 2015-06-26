#define MAX_SHORT 32767

//uniform hash function- Bob Jenkins' One-at-a-time algorithm
uint hash( uint x ) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}

//construct a float between 0 and 1 from the low 23 bits of a uint
float float01(uint m){
    const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    m |= ieeeOne;                          // Add fractional part to 1.0

    float  f = as_float( m );       // Range [1:2]
    return f - 1.0;                        // Range [0:1]
}

// Pseudo-random value in half-open range [0:1].
float random( float x ) { return float01(hash(as_uint(x))); }

//simple wrapper function
float randNext(float* randPrev){
	*randPrev = random(*randPrev);
	return (*randPrev);
}

//checks if a square of size s will fit at the current location
bool fitSquare(global char* puzzles, const char height, const char width, char row, char col, char size){
    //check the square is within the bounds
    if(row+size-1>=height || col+size-1>=width || row < 0 || col < 0){
        return false;
    }

    //check if any squares within the size are -1
    uint id = get_global_id(0);
    for(char i = 0; i < size; i++){
        for(char j = 0; j < size; j++){
            //check if this square is -1
            if(puzzles[id * height * width + (row+i) * width + (col+j)]== -1){
                return false;
            }
        }
    }

    return true;
}

//sets each puzzle square to the maximum size square possible at that point
char setMaxes(global char* puzzles, const char height, const char width){
    const uint id = get_global_id(0);
    char s;
    char largest = 0;

    //calculate the maximum possible size at each position in the puzzle
    for(char row = 0; row < height; row++){
        for(char col = 0; col < width; col++){
            //ignore blocked squares
            if(puzzles[id * height * width + row * width + col] != -1){
                //correctly initialise s
                if(puzzles[id * height * width + row * width + col] == 0){
                    s = 1;
                }else{
                    s = puzzles[id * height * width + row * width + col];
                }

                //keep incrementing s until we find a size of square that won't fit
                while(fitSquare(puzzles,height,width,row,col,s+1)==true){
                    s++;
                }

                //set the value in the puzzle and update the starting point for the surrounding squares
                for(char x = 0; x < s; x++){
                    for(char y = 0; y < s; y++){
                        puzzles[id * height * width + (row + x) * width + (col + y)] = s - max(x,y);
                    }
                }

                //check if this is the largest seen value
                if(s > largest){
                    largest = s;
                }
            }
        }
    }
    return largest;
}

//updates the surrounding maximum values after placing a square
void updateMaxes(global char* puzzles, const char height, const char width, char trow, char tcol, char s, short* squares){
    const uint id = get_global_id(0);

    //update the maximums of the surrounding squares that may be effected by this one
    //rows directly above:
    for(char row = max(trow-s+1,0); row < trow; row++){
        for(char col = tcol; col < tcol + s; col++){
            //check if this square is effected
            //if it is equal to s the square count must be updated
            if(puzzles[id * height * width + row * width + col] == s){
                (*squares)--;
            }

            //if the square is above zero it is an unallocated maximum, so must be updated
            //if the square is less than -1 it has already been allocated and so is not relevant
            if(puzzles[id * height * width + row * width + col] > 0){
                //check if it overlaps the new square
                if(row+puzzles[id * height * width + row * width + col] > trow){
                    //if there is an overlap recalculate the max size
                    puzzles[id * height * width + row * width + col] = trow - row;
                }
            }
        }
    }

    //columns directly to the left
    for(char row = trow; row < trow + s; row++){
        for(char col = max(tcol-s+1,0); col < tcol; col++){
            //check if this square is effected
            //if it is equal to s the square count must be updated
            if(puzzles[id * height * width + row * width + col] == s){
                (*squares)--;
            }

            //if the square is above zero it is an unallocated maximum, so must be updated
            //if the square is less than -1 it has already been allocated and so is not relevant
            if(puzzles[id * height * width + row * width + col] > 0){
                //check if it overlaps the new square
                if(col+puzzles[id * height * width + row * width + col] > tcol){
                    //if there is an overlap recalculate the max size
                    puzzles[id * height * width + row * width + col] = tcol - col;
                }
            }
        }
    }
    //columns and rows diagonally above on the left
    for(char row = max(trow-s+1,0); row < trow; row++){
        for(char col = max(tcol-s+1,0); col < tcol; col++){
            //check if this square is effected
            //if it is equal to s the square count must be updated
            if(puzzles[id * height * width + row * width + col] == s){
                (*squares)--;
            }

            //if the square is above zero it is an unallocated maximum, so must be updated
            //if the square is less than -1 it has already been allocated and so is not relevant
            if(puzzles[id * height * width + row * width + col] > 0){
                //check if it overlaps the new square
                //we need to check rows and columns for the diagonal
                if(col+puzzles[id * height * width + row * width + col] > tcol && row+puzzles[id * height * width + row * width + col] > trow){
                    //if there is an overlap recalculate the max size
                    puzzles[id * height * width + row * width + col] = max(tcol - col,trow - row);
                }
            }
        }
    }
}

//complete a puzzle
void doPuzzle(global short* lengths, global char* puzzles, global char* solutions, const char height, const char width, short* len, float temp, float* randPrev){
    //global invocation number
    const uint id = get_global_id(0);

    //other variables
    char maxSquare;
    short squares;
    short squaresCount;
    short i;
    short count;
    char trow;
    char tcol;
    float prob;

    //get largest possible value for each square and largest value overall
    maxSquare = setMaxes(puzzles,height,width);

    //loop through from largest to size 2
    for(char s = maxSquare; s >= 2; s--){
        //first count the number of possible squares of this size
        squaresCount = 0;
        for(char row = 0; row < height; row++){
            for(char col = 0; col < width; col++){
                if(puzzles[id * height * width + row * width + col] == s){
                    squaresCount++;
                }
            }
        }

        //loop until all squares have been covered
        squares = squaresCount;
        while(squares > 0){
            //randomly select a square index
            i = convert_short(squares * randNext(randPrev));

            //find the square of size s with the index i
            count = 0;
            for(char row = 0; row < height && count<=i; row++){
                for(char col = 0; col < width && count<=i; col++){
                    if(puzzles[id * height * width + row * width + col] == s){
                        if(count==i){
                            trow = row;
                            tcol = col;
                        }
                        count++;
                    }
                }
            }

            //with a random probability do NOT fill the square with this size
            //this allows the square to be filled with smaller squares
            //prioritise larger rarer squares
            //also use simulated annealing so the process becomes more greedy over time
            //only perform this if s>2, as a square of size 2 will only filled with size 1
            prob = convert_float(s-2) / convert_float(maxSquare-2) * (1 / convert_float(squaresCount+1)) * temp;
            if(s > 2 && randNext(randPrev) < prob){
                //leave this square for smaller squares
                //in order to do this the max square size as size-1
                puzzles[id * height * width + trow * width + tcol] = s-1;
                //update square count
                squares--;
            }else{
                //fill this square with s
                //update the surrounding maximums and square count
                updateMaxes(puzzles,height,width,trow,tcol,s,&squares);               

                //add this square to the solution
                for(char x = 0; x < s; x++){
                    for(char y = 0; y < s; y++){
                        //check if we are overwriting another potential square of this size and if so update count
                        if(puzzles[id * height * width + (trow+x) * width + (tcol+y)] == s){
                            squares--;
                        }

                        //set this square to the negative of s so it does not interfere with subsequent rounds
                        //(this way we can tell apart maximums and placed squares)
                        //we will have to make this positive this later (the things I do for memory)
                        puzzles[id * height * width + (trow+x) * width + (tcol+y)] = -s;
                    }
                }
                
                //record that we have added to the solution length
                (*len)++;
            }
        }
    }   

    //now fill in the squares of size 1
    //we also need to negate all filled squares (as we negated them to differentiate from maximums earlier)
    for(char row = 0; row < height; row++){
        for(char col = 0; col < width; col++){
            //negate filled squares
            if(puzzles[id * height * width + row * width + col] < -1){
                puzzles[id * height * width + row * width + col] = -puzzles[id * height * width + row * width + col];
            }

            //remaining squares are already size 1, time to count them
            if(puzzles[id * height * width + row * width + col] == 1){
                (*len)++;
            }
        }
    } 
}

//aggregate the results across the workgroup
void aggregateWorkGroup(global short* lengths, global short* groupLengths,  global char* puzzles, global char* solutions, const char height, const char width){
    const uint id = get_global_id(0);

    //find best result in the group
    barrier(CLK_GLOBAL_MEM_FENCE);

    short bestLength = MAX_SHORT;
    uint bestIndex;
    //check the others solutions in the work group
    for(uint gid = get_group_id(0) * get_local_size(0); gid < get_group_id(0) * get_local_size(0) + get_local_size(0); gid++){
        //see if this invocation has a better result
        if(lengths[gid] < bestLength){
            bestIndex = gid;
            bestLength = lengths[gid];
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    //check if this is the best result yet
    if(bestLength < groupLengths[get_group_id(0)]){
        //the invocation with the best result copies it to the solution
        //the other invocations copy the best result from the invocation that calculated it
        if(bestIndex == id){
            //update best length
            groupLengths[get_group_id(0)] = bestLength;

            //update best solution
            for(char row = 0; row < height; row++){
                for(char col = 0; col < width; col++){
                    solutions[get_group_id(0) * height * width + row * width + col] = puzzles[bestIndex * height * width + row * width + col];
                }
            }
        }else{
            //otherwise copy the solution from the best index
            for(char row = 0; row < height; row++){
                for(char col = 0; col < width; col++){
                    puzzles[id * height * width + row * width + col] = puzzles[bestIndex * height * width + row * width + col];
                }
            }
        }
    }else{
        //if this is not the best result, all invocations copy the best result from the solution
        for(char row = 0; row < height; row++){
            for(char col = 0; col < width; col++){
                puzzles[id * height * width + row * width + col] = solutions[get_group_id(0) * height * width + row * width + col];
            }
        }
    }

    //ensure all invocations have the best puzzle
    barrier(CLK_GLOBAL_MEM_FENCE);
}

void perturbPuzzle(global char* puzzles, global char* solutions, const char height, const char width, short* len, float* randPrev){
    //randomly set some squares to -1 in the puzzle
    //these squares will be filled in with the solution from the best solution later
    //also count the number of squares we set to zero and put in the lengths buffer
    const uint id = get_global_id(0);
   
    char s;
    (*len) = 0;
    float prob = 0.5;
    for(char row = 0; row < height; row++){
        for(char col = 0; col < width; col++){
            s = puzzles[id * height * width + row * width + col];
            //ignore blocked squares
            if(s != -1){
                //always set size 1 squares to 0 not -1 as they are trivial to fill
                if(randNext(randPrev) > prob && s > 1){
                    //set entire square to -1
                    (*len)++;
                    for(int x = 0; x < s; x++){
                        for(int y = 0; y < s; y++){
                            puzzles[id * height * width + (row+x) * width + (col+y)] = -1;
                        }
                    }
                }else{
                    //set entire square to 0
                    for(int x = 0; x < s; x++){
                        for(int y = 0; y < s; y++){
                            puzzles[id * height * width + (row+x) * width + (col+y)] = 0;
                        }
                    }
                }
            }
        }
    }  
}

//recover a perturbed puzzle to a full puzzle
void recoverPuzzle(global char* puzzles, global char* solutions, const char height, const char width){
    const uint id = get_global_id(0);

    char s;
    for(char row = 0; row < height; row++){
        for(char col = 0; col < width; col++){
            //if the square is -1 in the puzzle but not -1 in the solution, it was randomly set to -1 as a peturbation
            if(puzzles[id * height * width + row * width + col] == -1 && solutions[get_group_id(0) * height * width + row * width + col] != -1){
                //therefore we need to copy this square from the solution
                s = solutions[get_group_id(0) * height * width + row * width + col];

                for(char x = 0; x < s; x++){
                    for(char y = 0; y < s; y++){
                        puzzles[id * height * width + (row+x) * width + (col+y)] = s;
                    }
                }
            }
        }
    }
}

//main
__kernel void montecarlo(
	__global short* lengths, 
	__global short* groupLengths, 
	__global char* puzzles, 
	__global char* solutions,
	const char height, 
	const char width, 
	const int iterationNo){

    //global invocation number
    const uint id = get_global_id(0);

    //simulated annealing coefficients
    const float alpha = 0.95;
    float temp = 1.0;

    //other iteration variables
    float randPrev;
    short len = 0;

    //first iteration is completely random

    //initialise random number generator with id and iteration number
    randPrev = float01(hash(as_uint(convert_float(id)+1) ^ hash(as_uint(convert_float(0)))));

    //do the puzzle
    doPuzzle(lengths,puzzles,solutions,height,width,&len,temp,&randPrev);

    //copy to lengths buffer
    lengths[id] = len;

    //aggregate the best result across the workgroup
    aggregateWorkGroup(lengths, groupLengths, puzzles, solutions, height, width);

    //now refine the puzzle by perturbing the best solution and solving for that for the remaining iterations
	for(int iter=1;iter<iterationNo;iter++){

        //seed random number generator using both iteration number and invocation number
        //this ensures future iterations will have different results
        randPrev = float01(hash(as_uint(convert_float(id)+1) ^ hash(as_uint(convert_float(iter)))));

        //simulated annealing temperature
        temp *= alpha;

        //perturb the puzzle
        perturbPuzzle(puzzles,solutions,height,width,&len,&randPrev);

	    //fill the puzzle with squares and get solution length
	    doPuzzle(lengths,puzzles,solutions,height,width,&len,temp,&randPrev);
	    
        //recover the puzzle (i.e. to its unperturbed form)
        recoverPuzzle(puzzles,solutions,height,width);

	    //copy to lengths buffer
	    lengths[id] = len;

        //aggregate the best result across the workgroup
        aggregateWorkGroup(lengths, groupLengths, puzzles, solutions, height, width);
	}
}