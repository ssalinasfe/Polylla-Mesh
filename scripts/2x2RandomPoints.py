# Generate random points inside a square 2x2 
# Input: nPoints: number of points inside the square and in the borders
# Output: A .node file with the random points

import random
import sys
import math

def write_points(nPoints, xPoints, yPoints):
    largo =len(xPoints)
    f = open('2x2_' + str(largo) + '.node', 'w')
    #f = open('autodata.node', 'w')
    f.write("{} 2 0 0\n".format(largo))
    for i in range(0,largo):
        f.write('{0} {1} {2}\n'.format(i, xPoints[i], yPoints[i]))
    f.write('\n')
    f.close()

def generate_random_points(nPoints, xPoints, yPoints, tolerance):
    uwu = 2.0
    s = set([(0, 0), (2,2), (0, 2), (2,0)])

    while len(s) != nPoints:
        x = random.uniform(0.000000000001, 2.00000000 )
        y = random.uniform(0.000000000001, 2.00000000 )
        #print("premove ",x,y)
        s.add(move_point(uwu, x, y, tolerance))
        #print("postmove ",x,y)
        #print(len(s))

    for _ in range(0,len(s)):
        l = s.pop()
        xPoints.append(l[0])
        yPoints.append(l[1])


def move_point(nPoints, xPoints, yPoints, tolerance):
    #tolerance =  0.0001  
    #tolerance =  5/(nPoints/10)
    n = nPoints
    r = random.uniform(0, 1)
    if r > 0.5:
        if xPoints >= nPoints*(1.0-tolerance): 
            xPoints = n
        if yPoints >= nPoints*(1.0-tolerance): 
            yPoints = n
        if xPoints <= nPoints*tolerance: 
            xPoints = 0
        if yPoints <= nPoints*tolerance: 
            yPoints = 0
    else:
        if xPoints <= nPoints*tolerance: 
            xPoints = 0            
        if yPoints <= nPoints*tolerance: 
            yPoints = 0
        if xPoints >= nPoints*(1.0-tolerance): 
            xPoints = n
        if yPoints >= nPoints*(1.0-tolerance): 
            yPoints = n
    #print("returning", xPoints, yPoints)
    return (xPoints, yPoints)

#Input Largo inicial, incremento

def main():
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]

    nPoints = int(argument_list[0])
    #tolerance = float(argument_list[1])
    tolerance = 1/10**(0.5*math.log10(nPoints) + 0.5)
    tolerance *= 2
    random.seed(138)
    xPoints = []
    yPoints = []
    
    generate_random_points(nPoints, xPoints, yPoints, tolerance)
    write_points(nPoints, xPoints, yPoints)
    

if __name__ == "__main__":
    main()
