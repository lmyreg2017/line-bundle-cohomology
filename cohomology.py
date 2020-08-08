#########################################################################
#                                                                       #
# Cohomology of line bundles on toric varieties using the Cech complex. #
#                                                                       #
#########################################################################

from subprocess import Popen, PIPE
import os
import tempfile
import numpy as np

def choose_k_in_s(k, s):
    if k == 0:
        yield []
        return
    
    for i in xrange(0,len(s)-(k-1)):
        for combination in choose_k_in_s(k-1, s[i+1:]):
            yield [s[i]]+combination


def compute_box(rays, divisor):
    dim = len(rays[0])
    # min, max
    boxsize = [(0,0)]*dim
	
    for ray_system in choose_k_in_s(dim, zip(rays, divisor)):
        rays = np.array( [ item[0] for item in ray_system ] )
        a = np.array( [ -item[1] for item in ray_system ] )
        #rays = matrix(CC, [item[0] for item in ray_system])
        #a = matrix(ZZ, [-item[1] for item in ray_system])
        print( rays )
        try:
            m = np.linalg.inv( rays )
        except:
        #except ZeroDivisionError:
            # No solution, just try the next collection
            continue

        solution = m * a.transpose()

        for i in range(dim):
            #m_i = int(RR(real(solution[i][0])))
            m_i = int(solution[i][0])
            min_m, max_m = boxsize[i]
            min_m = min(min_m, m_i)
            max_m = max(max_m, m_i)
            boxsize[i] = (min_m, max_m)
		
    for i in range(dim):
        boxsize[i] = (boxsize[i][0]-1,boxsize[i][1]+1)

    return boxsize

def compute_kth_cohomology(rays, cones, divisor, k):
    # Some basic sanity checks.
    assert len(divisor) == len(rays)
    assert all(all(ray < len(rays) for ray in cone) for cone in cones)
    assert len(cones) > len(rays[0])

    # Dimension of the N lattice. We will compute up to C^{dim}, d^{dim}
    dim = len(rays[0])

    # A bit more convenient representation for computing intersections.
    cones = [set(cone) for cone in cones]

    def dot(a,b):
        return sum([i*j for i,j in zip(a,b)])

    # Pass the relevant information to the C code
    fname = open("demofile.txt", "w")
    #fname = tempfile.TemporaryFile()
    fd = fname
    #fd = open(fname, "wb")
    box = compute_box(rays, divisor)
    
    # Definition of the box.
    fd.write(str(len(box))+"\n")
    for interval in box:
        fd.write("%d %d\n" % interval)

    # The rays.
    fd.write(str(len(rays))+"\n")
    for ray in rays:
        for xi in ray:
            fd.write("%d " % (xi,))
        fd.write("\n")

    # The divisor.
    for ai in divisor:
        fd.write("%d " % (ai,))
    fd.write("\n")

    # The cones.
    fd.write(str(len(cones))+"\n")
    for cone in cones:
        fd.write("%d\n" % (len(cone),))
        for ray in cone:
            fd.write("%d " % (ray,))
        fd.write("\n")
    
    fd.close()

    print( str( k ) )
    print( fname )
    cech = Popen(['cech_cohomology', fname, str(k)], stdout=PIPE)
    
    output, errors = cech.communicate()
    result = int(output.strip())

    os.unlink(fname)

    return result






######################################################################
# End of the code for computing cohomology.                          #
# The following lines are examples.                                  #
######################################################################

## dP_1
#rays = [[1,0], [0,1], [-1,-1], [0,-1]]
#cones = [[0,1], [1,2], [2,3], [3,0]]
## Divisor we are interested in, in terms of the divisors associated
## with the rays in the fan. We take for example D = 5D_1 - 2D_4
#D = [5, 0, 0, -2]

## P^2
#rays = [[1,0],[0,1],[-1,-1]]
#cones = [[0,1], [1,2], [2,0]]
#D = [7,1,-1]

## An example from 1002.1894 (around eq. 48)
#rays=[[0,0,0,1,0],[0,0,0,0,1],[0,0,0,-2,-3],[-1,-1,-1,-8,-12],[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0]]
#cones=[[0,1,3,4,5],[0,1,3,4,6],[0,1,3,5,6],[0,1,4,5,6],[0,2,3,4,5],[0,2,3,4,6],[0,2,3,5,6],[0,2,4,5,6],[1,2,3,4,5],[1,2,3,4,6],[1,2,3,5,6],[1,2,4,5,6]]
#D=[-1,-1,-1,-1,-1,-1,-5]

#for k in range(0,len(rays[0])+1):
    #print "Computing H^{%d}" % (k,)
    #Hk = compute_kth_cohomology(rays, cones, D, k)
    #print "H^{%d} = %d" % (k,Hk)

def add(x,y):
	return [a+b for a,b in zip(x,y)]

def minus(x):
	return [-1*i for i in x]

rays=[(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0), (0, 0, 0,1, 0, 0), (0, 0, 0, 0, 1, 0), (-3, -2, 0, 0, 0, 0), (6, 4, 1, 1, 1, 0),(-6, -4, 0, -1, 0, 0), (0, 0, 2, 1, 1, 3), (-3, -2, -2, -1, -1, -2), (9,6, 2, 1, 1, 0)]
conesbase=[[1,2,3,4],[1,2,3,8],[1,2,4,7],[1,2,5,7],[1,2,5,8],[1,3,4,6],[1,3,6,8],[1,4,6,7],[1,5,6,7],[1,5,6,8],[2,3,4,6],[2,3,6,8],[2,4,6,7],[2,5,6,7],[2,5,6,8]]
cones=[]
for cone in conesbase:
         for extension in ([9,10],[9,11],[10,11]):
                 cones.append([n-1 for n in cone+extension])
NstarBase = [-1,-1,-1,-1,0,-1,0,-1,0,0,0] # = -6I
NstarBaseShort = [-2,0,0,0,0,0,0,0,0,0,0]
NstarDiv = [0,0,0,0,-1,0,0,0,0,0,0]
NstarFiber = [0,0,0,0,0,0,0,0,0,-2,0]

# Next  3 are to compute H^2(N*)
#D = NstarBase
#print "H^2 NstartBase"
#print compute_kth_cohomology(rays,cones,D,2)

#D = NstarFiber
#print "H^2 NstarFiber"
#print compute_kth_cohomology(rays,cones,D,2)

#D = NstarDiv
#print "H^2 NstartDiv"
#print compute_kth_cohomology(rays,cones,D,2)

##Computing H^4(wedge^2N*) via Serre duality

#D = add(add(add([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],minus(NstarBaseShort)),minus(NstarFiber)),minus(NstarDiv))
#print D
#print "H^2(K tensor wedge^3N*)"
#print compute_kth_cohomology(rays,cones,D,2)

## Canonical line bundle
canon=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]

D=add(NstarBase,NstarFiber)
D1=add(NstarDiv,NstarFiber)
D2=add(NstarBase,NstarDiv)

print( D )
print( D1 )
print( D2 )

## The third cohomology is the slowest to compute. The following code
## takes a while to finish.
print( compute_kth_cohomology(rays,cones,D,3) )

print( compute_kth_cohomology(rays,cones,D1,3) )

print( compute_kth_cohomology(rays,cones,D2,3) )
