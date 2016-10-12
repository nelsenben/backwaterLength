__author__ = "Ben Nelsen"

import numpy as np
import os

class Channel:

    def __init__(self):
        self.flow = 0
        self.channelwidth = 0
        self.mannings = 0
        self.epsilon = 0
        self.S01 = 0
        self.y01 = 0
        self.S2 = 0
        self.Fr1 = 0
        self.S02 = 0
        self.y02 = 0
        self.S1 = 0
        self.Fr2 = 0
        self.jump = False
        self.jumpLocation = 'None'
        self.depthi = 0
        self.depthf = 0
        self.slope = 0


    def readConfigurationfile(self, s_folderPath, s_fileName):

        s_filePath = os.path.join(s_folderPath, s_fileName)
        f_data = file(s_filePath, 'r')

        for line in f_data:
            column = line.split('\t')
            if column[0] == "Q":
                self.flow = float(column[1])
            elif column[0] == "b":
                self.channelwidth = float(column[1])
            elif column[0] == "n":
                self.mannings = float(column[1])
            elif column[0] == "eps":
                self.epsilon = float(column[1])
            elif column[0] == "S01":
                self.S01 = float(column[1])
            elif column[0] == "S02":
                self.S02 = float(column[1])
            else:
                print('value not understood...')

    def manningsEquation(self, slope, mannings, width, normalDepth):

        flow = 1.49/mannings*(width*normalDepth)*(width*normalDepth/(width+2*normalDepth))**(2/3)*slope**0.5
        return flow

    def sequentDepths(self, normalDepth, unitFlow):

        depthSequent = normalDepth/2*(np.sqrt(1+8*unitFlow**2/(32.2*normalDepth**3))-1)
        return depthSequent

    def calculateFroude(self, velocity, depth):

        froudeNumber = velocity/np.sqrt(depth*32.2)
        return froudeNumber

    def lengthItteration(self):
        d_totalL1 = 0
        d_totalL2 = 0
        changeTotal1 = 9999
        itteration = 1

        while not (changeTotal1 < self.epsilon*d_totalL2):
            depth1 = self.depthi
            velocity1 = self.flow / (self.depthi * self.channelwidth)
            energy1 = self.depthi + velocity1 ** 2 / (2 * 32.2)
            radius1 = self.channelwidth * self.depthi / (self.channelwidth + 2 * self.depthi)

            d_totalL1 = d_totalL2
            d_totalL2 = 0
            unitY = (self.depthf-self.depthi)/itteration
            for index in range(0,itteration):
                y2 = depth1 + unitY
                velocity2 = self.flow/(y2*self.channelwidth)
                energy2 = y2 + velocity2 ** 2 / (2 * 32.2)
                radius2 = self.channelwidth * y2 / (self.channelwidth + 2 * y2)
                averageRadius = (radius1 + radius2)/2
                averageVelocity = (velocity1 + velocity2)/2
                averageSlope = (averageVelocity*self.mannings/(1.49*averageRadius**(2/3)))**2
                changeLength = (energy2-energy1)/(self.slope-averageSlope)
                d_totalL2 = d_totalL2 + changeLength

                depth1 = y2
                velocity1 = velocity2
                energy1 = energy2
                radius1 = radius2

            changeTotal1 = abs(d_totalL2-d_totalL1)

        print('debug')






    def backwaterDetermination(self):
        for index in range(0,2):
            if index == 0:
                slope = self.S01
            elif index == 1:
                slope = self.S02
            else:
                print('error')
            depthTry = 1
            flowTry = 0
            depthBetter = 1
            while not (abs(self.flow-flowTry) < self.epsilon*self.flow):
                depthTry = depthBetter
                flowTry = self.manningsEquation(slope, self.mannings, self.channelwidth, depthTry)
                depthBetter = depthTry*(self.flow/flowTry+1)/2
            depthSequent = self.sequentDepths(depthTry, self.flow/self.channelwidth)
            froudeNumber = self.calculateFroude(self.flow/(self.channelwidth*depthTry),depthTry)
            if index == 0:
                self.y01 = depthTry
                self.S2 = depthSequent
                self.Fr1 = froudeNumber

            elif index == 1:
                self.y02 = depthTry
                self.S1 = depthSequent
                self.Fr2 = froudeNumber

            else:
                print('error')

        if (self.Fr1 > 1 or self.Fr2 > 1) and (self.Fr1 < 1 or self.Fr2 < 1):
            self.jump = True

        if self.jump:
            if self.S2 > self.y02:
                self.jumpLocation = "downstream"
                self.depthi = self.y01
                self.depthf = self.S1
                self.slope = self.S02

            elif self.S2 < self.y02:
                self.jumpLocation = "upstream"
                self.depthi = self.S2
                self.depthf = self.y02
                self.slope = self.S01

            else:
                print('error')
        else:
            print('No Jump occurs...')

        self.lengthItteration()

