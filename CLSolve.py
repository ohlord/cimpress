import pyopencl as cl
import numpy as np
import os
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

class CLSolver:
    
    def loadKernel(self, filename):
        #read in the OpenCL source file as a string
        f = open(filename, 'r')
        contents = "".join(f.readlines())
        
        #create the program
        self.program = cl.Program(self.context, contents).build()
        
    def init(self, filename="montecarlo.cl"):  
        #setup context and command queue
        self.platform = cl.get_platforms()
        self.gpu_devices = self.platform[0].get_devices(device_type = cl.device_type.GPU)
        
        #print device details
        for i in range(len(self.gpu_devices)):
            device = self.gpu_devices[i]
            print "Device:", i
            print device
            print "Global Memory:", device.get_info(cl.device_info.GLOBAL_MEM_SIZE) / 1000, "KB"
            print "Local Memory:",device.get_info(cl.device_info.LOCAL_MEM_SIZE) / 1000, "KB"
        
        #create openCL context and queue
        self.context = cl.Context(devices=self.gpu_devices)
        self.queue = cl.CommandQueue(self.context)
        
        #load kernel from file
        self.loadKernel("Solvers/montecarlo.cl")
    
    def initBuffers(self,puzzle):
        #define lengths buffer and copy to the GPU
        #as we will not read from this buffer later, mapping is not required
        self.lengths = np.full(self.simulations,np.iinfo(np.int16).max,dtype=np.int16)
        self.lengthsBuffer = cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.lengths)
         
        #define buffer for aggregated lengths for each workgroup
        self.groupLengths = np.full(self.workGroups,np.iinfo(np.int16).max,dtype=np.int16)
        self.groupLengthsBuffer = cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.USE_HOST_PTR, hostbuf=self.groupLengths)
        
        #map group lengths buffer
        cl.enqueue_map_buffer(self.queue,self.groupLengthsBuffer,cl.map_flags.READ,0,self.groupLengths.shape,self.groupLengths.dtype)
        
        #get the input puzzle ready for the kernel; convert to 8 bit int (char)
        p = np.array(puzzle['puzzle']).astype(np.int8)
        #subtract 1 so that -1 denotes a gap and 0 denotes a square to be filled
        p = p - np.ones_like(p,dtype=p.dtype)
        
        #copy the puzzle, one for each simulation
        self.puzzles = np.zeros((self.simulations,self.height,self.width),dtype=p.dtype)
        self.puzzles[:,0:self.height,0:self.width] = p
    
        #define puzzles buffer and copy data (we do not need to worry about getting data out of this buffer, so mapping isn't required)
        #this buffer contains the input puzzles, one for each invocation (the puzzle is too large to hold in local or shared memory)
        self.puzzlesFlattened = self.puzzles.ravel()
        self.puzzlesBuffer = cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.puzzlesFlattened)
        
        #define output buffer for best solutions aggregated across workgroups
        self.solutions = self.puzzles[0:self.workGroups]
        self.solutionsFlattened = self.solutions.ravel()
        self.solutionsBuffer = cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.USE_HOST_PTR, hostbuf=self.solutionsFlattened)

        #map solutions buffer
        cl.enqueue_map_buffer(self.queue,self.solutionsBuffer,cl.map_flags.READ,0,self.solutionsFlattened.shape,self.solutions.dtype)
        
    
    def solve(self,puzzle,simulations = 16384, iterations = 35, workGroupSize = 128):
        self.simulations = simulations
        self.iterations = iterations
        self.workGroupSize = workGroupSize
        self.workGroups = int(self.simulations / self.workGroupSize)
        self.width = np.int8(puzzle['width'])
        self.height = np.int8(puzzle['height'])
        
        #initialise buffers
        self.initBuffers(puzzle)
        
        #create kernel
        self.kernel = cl.Kernel(self.program,"montecarlo")
        self.kernel.set_args(self.lengthsBuffer,self.groupLengthsBuffer,self.puzzlesBuffer,self.solutionsBuffer,self.height,self.width,np.int32(self.iterations))
        
        #execute program for a number of iterations
        cl.enqueue_nd_range_kernel(self.queue,self.kernel,(self.simulations,),(self.workGroupSize,))
        
        #unmap group lengths buffer from device
        cl.enqueue_map_buffer(self.queue,self.groupLengthsBuffer,cl.map_flags.WRITE,0,self.groupLengths.shape,self.groupLengths.dtype)
        self.groupLengths = self.groupLengthsBuffer.get_host_array(self.groupLengths.shape,dtype=self.groupLengths.dtype)

        #unmap solutions buffer from device
        cl.enqueue_map_buffer(self.queue,self.solutionsBuffer,cl.map_flags.WRITE,0,self.solutionsFlattened.shape,self.solutions.dtype)
        self.solutions = self.solutionsBuffer.get_host_array(self.solutions.shape,dtype=self.solutions.dtype)
        
        #release buffers
        self.lengthsBuffer.release()
        self.groupLengthsBuffer.release()
        self.puzzlesBuffer.release()
        self.solutionsBuffer.release()

        #get the best solution
        i = self.groupLengths.argmin()
        bestSolution = np.array(self.solutions[i])
        
        #convert solution to list format used by challenge
        solution = []
        for row in range(0,puzzle['height']):
            for col in range(0,puzzle['width']):
                if bestSolution[row][col]!=-1:
                    s = bestSolution[row][col]
                    
                    #add to solution list
                    solution.append({'X': int(col),'Y': int(row),'Size':int(s)})
                    
                    #clear cells in solution
                    for i in range(0,s):
                        for j in range(0,s):
                            bestSolution[row+i][col+j]=-1
        
        return solution

