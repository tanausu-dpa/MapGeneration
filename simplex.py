# -*- coding: utf-8 -*-

######################################################################
######################################################################
######################################################################
#                                                                    #
# simplex.py                                                         #
#                                                                    #
# Tanaus\'u del Pino Alem\'an                                        #
#   Instituto de Astrof\'isica de Canarias                           #
#                                                                    #
######################################################################
######################################################################
#                                                                    #
# Class with simplex noise generation                                #
#                                                                    #
######################################################################
######################################################################
#                                                                    #
# This program is free software: you can redistribute it and/or      #
# modify it under the terms of the GNU General Public License as     #
# published by the Free Software Foundation, either version 3 of     #
# the License, or (at your option) any later version.                #
#                                                                    #
# This program is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                           #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with this program.  If not, see                              #
# <https://www.gnu.org/licenses/>.                                   #
#                                                                    #
######################################################################
######################################################################
#                                                                    #
#  02/11/2018 - V1.0.0 - Added license. (TdPA)                       #
#                                                                    #
#  16/04/2018 - V0.0.0 - Start Code. (TdPA)                          #
#                                                                    #
######################################################################
######################################################################
######################################################################
#                                                                    #
# This code is direct translation of the code by Eliot Eshelman,     #
# with the following Copyright:                                      #
#                                                                    #
#* Copyright (c) 2007-2012 Eliot Eshelman                            #
#*                                                                   #
#* This program is free software: you can redistribute it and/or     #
#* modify it under the terms of the GNU General Public License as    #
#* published by the Free Software Foundation, either version 3 of    #
#* the License, or at your option) any later version.                #
#*                                                                   #
#* This program is distributed in the hope that it will be useful,   #
#* but WITHOUT ANY WARRANTY; without even the implied warranty of    #
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      #
#* GNU General Public License for more details.                      #
#*                                                                   #
#* You should have received a copy of the GNU General Public License #
#* along with this program. If not, see                              #
#* <http://www.gnu.org/licenses/>.                                   #
#*                                                                   #
# And the following comments:                                        #
#                                                                    #
#* 2D, 3D and 4D Simplex Noise functions return 'random' values in   #
#* (-1, 1).                                                          #
#*                                                                   #
#* This algorithm was originally designed by Ken Perlin, but my code #
#* has been adapted from the implementation written by Stefan        #
#* Gustavson (stegu@itn.liu.se)                                      #
#*                                                                   #
#* Raw Simplex noise functions return the value generated by Ken's   #
#* algorithm.                                                        #
#*                                                                   #
#* Scaled Raw Simplex noise functions adjust the range of values     #
#* returned from the traditional (-1, 1) to whichever bounds are     #
#* passed to the function.                                           #
#*                                                                   #
#* Multi-Octave Simplex noise functions compine multiple noise       #
#* values to create a more complex result. Each successive layer of  #
#* noise is adjusted and scaled.                                     #
#*                                                                   #
#* Scaled Multi-Octave Simplex noise functions scale the values      #
#* returned from the traditional (-1,1) range to whichever range is  #
#* passed to the function.                                           #
#*                                                                   #
#* In many cases, you may think you only need a 1D noise function,   #
#* but in practice 2D is almost always better.  For instance, if     #
#* you're using the current frame number as the parameter for the    #
#* noise, all objects will end up with the same noise value  at each #
#* frame. By adding a second parameter on the second dimension, you  #
#* can ensure that each gets a unique noise value and they don't all #
#* look identical.                                                   #
#*                                                                   #
#                                                                    #
######################################################################
######################################################################
######################################################################

import math

######################################################################
######################################################################
######################################################################
######################################################################

class simplex_class():
    ''' Class with functions to generate simplex/Perlin noise
    '''

######################################################################
######################################################################

    def __init__(self, octaves=1, persistence=1.0, scale=1.0):
        ''' Initialize class
        '''

        # Set defaults
        self.__octaves = octaves
        self.__persistence = persistence
        self.__scale = scale
        self.__init_grad()

######################################################################
######################################################################

    def octave_noise_35d(self,x,y,z,w,octaves=None,persistence=None, \
                         scale=None):
        ''' * 3.5D Multi-octave Simplex noise.
            * For each octave, a higher frequency/lower amplitude
            * function will be added to the original. 
            * The higher the persistence [0-1], the more of each
            * succeeding octave will be added.
        '''

        if octaves is None: octaves = self.__octaves
        if persistence is None: persistence = self.__persistence
        if scale is None: scale = self.__scale

        total = 0.0
        amplitude = 1.0

        try:
            frequency = float(scale)
        except ValueError:
            return 1
        except:
            return -1

        #* We have to keep track of the largest possible amplitude,
        #* because each octave adds more, and we need a value in
        #* [-1, 1].
        maxAmplitude = 0.0

        for ii in range(octaves):

            total += self.__raw_noise_4d(x*frequency, \
                                         y*frequency, \
                                         z*frequency, w)*amplitude
    
            frequency *= 2.0
            maxAmplitude += amplitude
            amplitude *= persistence

        return total/maxAmplitude
    
######################################################################
######################################################################

    def scaled_octave_noise_35d(self,x,y,z,w,octaves=None, \
                                     persistence=None, \
                                     scale=None, \
                                     loBound=0.0, \
                                     hiBound=1.0):
        ''' * 3.5D Scaled Multi-octave Simplex noise.
            * Returned value will be between loBound and hiBound.
        '''

        if octaves is None: octaves = self.__octaves
        if persistence is None: persistence = self.__persistence
        if scale is None: scale = self.__scale

        return self.octave_noise_35d(x,y,z,w,octaves=octaves, \
                                           persistence=persistence, \
                                           scale=scale)* \
                                            (hiBound - loBound)/2 + \
                                            (hiBound + loBound)/2

######################################################################
######################################################################

    def __raw_noise_4d(self,x,y,z,w):
        ''' * 4D raw Simplex noise
        '''

        #* The skewing and unskewing factors are hairy again for the
        #* 4D case
        F4 = (math.sqrt(5.0)-1.0)/4.0;
        G4 = (5.0-math.sqrt(5.0))/20.0;

        #* Skew the (x,y,z,w) space to determine which cell of 24
        #* simplices we're in
        s = (x + y + z + w)*F4; #* Factor for 4D skewing
        i = self.__fastfloor(x + s)
        j = self.__fastfloor(y + s)
        k = self.__fastfloor(z + s)
        l = self.__fastfloor(w + s)
        t = (i + j + k + l)*G4 #* Factor for 4D unskewing
        X0 = i - t #* Unskew the cell origin back to (x,y,z,w) space
        Y0 = j - t
        Z0 = k - t
        W0 = l - t

        x0 = x - X0; #* The x,y,z,w distances from the cell origin
        y0 = y - Y0;
        z0 = z - Z0;
        w0 = w - W0;

        #* For the 4D case, the simplex is a 4D shape I won't even
        #* try to describe.
        #* To find out which of the 24 possible simplices we're in,
        #* we need to determine the magnitude ordering of x0, y0, z0
        #* and w0.
        #* The method below is a good way of finding the ordering of
        #* x,y,z,w and then find the correct traversal order for the
        #* simplex we're in.
        #* First, six pair-wise comparisons are performed between
        #* each possible pair of the four coordinates, and the results
        #* are used to add up binary bits for an integer index.
        c1 = 32 if (x0 > y0) else 0
        c2 = 16 if (x0 > z0) else 0
        c3 =  8 if (y0 > z0) else 0
        c4 =  4 if (x0 > w0) else 0
        c5 =  2 if (y0 > w0) else 0
        c6 =  1 if (z0 > w0) else 0
        c = c1 + c2 + c3 + c4 + c5 + c6

        #* simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in
        #* some order.
        #* Many values of c will never occur, since e.g. x>y>z>w makes
        #* x<z, y<w and x<w impossible. Only the 24 indices which have
        #* non-zero entries make any sense.
        #* We use a thresholding to set the coordinates in turn from
        #* the largest magnitude.
        #* The number 3 in the "simplex" array is at the position of
        #* the largest coordinate.
        i1 = 1 if (self.__simplex[c][0]>=3) else 0
        j1 = 1 if (self.__simplex[c][1]>=3) else 0
        k1 = 1 if (self.__simplex[c][2]>=3) else 0
        l1 = 1 if (self.__simplex[c][3]>=3) else 0

        #* The number 2 in the "simplex" array is at the second
        #* largest coordinate.
        i2 = 1 if (self.__simplex[c][0]>=2) else 0
        j2 = 1 if (self.__simplex[c][1]>=2) else 0
        k2 = 1 if (self.__simplex[c][2]>=2) else 0
        l2 = 1 if (self.__simplex[c][3]>=2) else 0

        #* The number 1 in the "simplex" array is at the second
        #* smallest coordinate.
        i3 = 1 if (self.__simplex[c][0]>=1) else 0
        j3 = 1 if (self.__simplex[c][1]>=1) else 0
        k3 = 1 if (self.__simplex[c][2]>=1) else 0
        l3 = 1 if (self.__simplex[c][3]>=1) else 0

        #* The fifth corner has all coordinate offsets = 1, so no need
        #* to look that up.

        #* Offsets for second corner in (x,y,z,w) coords
        x1 = x0 - i1 + G4
        y1 = y0 - j1 + G4
        z1 = z0 - k1 + G4
        w1 = w0 - l1 + G4

        #* Offsets for third corner in (x,y,z,w) coords
        x2 = x0 - i2 + 2.0*G4
        y2 = y0 - j2 + 2.0*G4
        z2 = z0 - k2 + 2.0*G4
        w2 = w0 - l2 + 2.0*G4

        #* Offsets for fourth corner in (x,y,z,w) coords
        x3 = x0 - i3 + 3.0*G4
        y3 = y0 - j3 + 3.0*G4
        z3 = z0 - k3 + 3.0*G4
        w3 = w0 - l3 + 3.0*G4

        #* Offsets for last corner in (x,y,z,w) coords
        x4 = x0 - 1.0 + 4.0*G4
        y4 = y0 - 1.0 + 4.0*G4
        z4 = z0 - 1.0 + 4.0*G4
        w4 = w0 - 1.0 + 4.0*G4

        #* Work out the hashed gradient indices of the five simplex
        #* corners
        ii = i & 255
        jj = j & 255
        kk = k & 255
        ll = l & 255

        gi0 = self.__perm[ii + self.__perm[jj + self.__perm[kk + \
                          self.__perm[ll]]]] % 32
        gi1 = self.__perm[ii + i1 + self.__perm[jj + j1 + \
                          self.__perm[kk + k1 + \
                          self.__perm[ll+l1]]]] % 32
        gi2 = self.__perm[ii + i2 + self.__perm[jj + j2 + \
                          self.__perm[kk + k2 + \
                          self.__perm[ll+l2]]]] % 32
        gi3 = self.__perm[ii + i3 + self.__perm[jj + j3 + \
                          self.__perm[kk + k3 + \
                          self.__perm[ll+l3]]]] % 32
        gi4 = self.__perm[ii + 1 + self.__perm[jj + 1 + \
                          self.__perm[kk + 1 + \
                          self.__perm[ll+1]]]] % 32

        #* Calculate the contribution from the five corners
        t0 = 0.6 - x0*x0 - y0*y0 - z0*z0 - w0*w0
        if(t0<0):
            n0 = 0.0
        else:
            t0 *= t0
            n0 = t0*t0*self.__dot(self.__grad4[gi0], x0, y0, z0, w0)

        t1 = 0.6 - x1*x1 - y1*y1 - z1*z1 - w1*w1
        if(t1<0):
            n1 = 0.0
        else:
            t1 *= t1
            n1 = t1*t1*self.__dot(self.__grad4[gi1], x1, y1, z1, w1)

        t2 = 0.6 - x2*x2 - y2*y2 - z2*z2 - w2*w2
        if(t2<0):
            n2 = 0.0
        else:
            t2 *= t2
            n2 = t2*t2*self.__dot(self.__grad4[gi2], x2, y2, z2, w2)

        t3 = 0.6 - x3*x3 - y3*y3 - z3*z3 - w3*w3
        if(t3<0):
            n3 = 0.0
        else:
            t3 *= t3
            n3 = t3*t3*self.__dot(self.__grad4[gi3], x3, y3, z3, w3)

        t4 = 0.6 - x4*x4 - y4*y4 - z4*z4 - w4*w4
        if(t4<0):
            n4 = 0.0;
        else:
            t4 *= t4
            n4 = t4*t4*self.__dot(self.__grad4[gi4], x4, y4, z4, w4)

        #* Sum up and scale the result to cover the range [-1,1]
        return 27.0*(n0 + n1 + n2 + n3 + n4)

######################################################################
######################################################################

    def __fastfloor(self,x):
        '''
        '''

        return int(x) if (x>0) else int(x) - 1

######################################################################
######################################################################

    def __dot(self,g,x,y,z,w):
        '''
        '''

        return g[0]*x + g[1]*y + g[2]*z + g[3]*w

######################################################################
######################################################################

    def __init_grad(self):
        ''' Initializes gradient variables
        '''

        #* The gradients are the midpoints of the vertices of a
        #* hypercube.
        self.__grad4 = [[0,1,1,1],[0,1,1,-1],[0,1,-1,1],[0,1,-1,-1], \
                    [0,-1,1,1],[0,-1,1,-1],[0,-1,-1,1],[0,-1,-1,-1], \
                        [1,0,1,1],[1,0,1,-1],[1,0,-1,1],[1,0,-1,-1], \
                    [-1,0,1,1],[-1,0,1,-1],[-1,0,-1,1],[-1,0,-1,-1], \
                        [1,1,0,1],[1,1,0,-1],[1,-1,0,1],[1,-1,0,-1], \
                    [-1,1,0,1],[-1,1,0,-1],[-1,-1,0,1],[-1,-1,0,-1], \
                        [1,1,1,0],[1,1,-1,0],[1,-1,1,0],[1,-1,-1,0], \
                    [-1,1,1,0],[-1,1,-1,0],[-1,-1,1,0],[-1,-1,-1,0]]


        #* Permutation table.  The same list is repeated twice.
        self.__perm = [151,160,137,91,90,15,131,13,201,95,96, \
                       53,194,233,7,225,140,36,103,30,69,142, \
                       8,99,37,240,21,10,23,190,6,148,247, \
                       120,234,75,0,26,197,62,94,252,219,203, \
                       117,35,11,32,57,177,33,88,237,149,56, \
                       87,174,20,125,136,171,168,68,175,74,165, \
                       71,134,139,48,27,166,77,146,158,231,83, \
                       111,229,122,60,211,133,230,220,105,92, \
                       41,55,46,245,40,244,102,143,54,65,25, \
                       63,161,1,216,80,73,209,76,132,187,208, \
                       89,18,169,200,196,135,130,116,188,159, \
                       86,164,100,109,198,173,186,3,64,52,217, \
                       226,250,124,123,5,202,38,147,118,126, \
                       255,82,85,212,207,206,59,227,47,16,58, \
                       17,182,189,28,42,223,183,170,213,119, \
                       248,152,2,44,154,163,70,221,153,101, \
                       155,167,43,172,9,129,22,39,253,19,98, \
                       108,110,79,113,224,232,178,185,112, \
                       104,218,246,97,228,251,34,242,193,238, \
                       210,144,12,191,179,162,241,81,51,145, \
                       235,249,14,239,107,49,192,214,31,181, \
                       199,106,157,184,84,204,176,115,121,50, \
                       45,127,4,150,254,138,236,205,93,222, \
                       114,67,29,24,72,243,141,128,195,78,66, \
                       215,61,156,180]*2

        #* A lookup table to traverse the simplex around a given point
        #* in 4D.
        self.__simplex = [[0,1,2,3],[0,1,3,2],[0,0,0,0], \
                          [0,2,3,1],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[1,2,3,0],[0,2,1,3], \
                          [0,0,0,0],[0,3,1,2],[0,3,2,1], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [1,3,2,0],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [1,2,0,3],[0,0,0,0],[1,3,0,2], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [2,3,0,1],[2,3,1,0],[1,0,2,3], \
                          [1,0,3,2],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[2,0,3,1],[0,0,0,0], \
                          [2,1,3,0],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [2,0,1,3],[0,0,0,0],[0,0,0,0], \
                          [0,0,0,0],[3,0,1,2],[3,0,2,1], \
                          [0,0,0,0],[3,1,2,0],[2,1,0,3], \
                          [0,0,0,0],[0,0,0,0],[0,0,0,0], \
                          [3,1,0,2],[0,0,0,0],[3,2,0,1], \
                          [3,2,1,0]]

######################################################################
######################################################################
######################################################################
######################################################################
