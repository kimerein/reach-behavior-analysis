# for all .AVI files in directory, convert to .MP4

import numpy as np
import cv2
import os
import sys

def convertLowSpeedAVIs():
    i=0
    useTheseFiles = []
    dirname = sys.argv[1]
    for file in os.listdir(dirname):
        if file.endswith(".AVI"):
            print(os.path.join(dirname, file))
            useTheseFiles.append(os.path.join(dirname, file))
            i=i+1

    for x in range(i):
        currFile = useTheseFiles[x]
        print("Now processing: " + currFile)
        
        cap = cv2.VideoCapture(currFile)

        maxFrames = 100
        frameCount = 0

        while(cap.isOpened()):
            ret, frame = cap.read()
            if frameCount%10 == 0:
                print(frameCount)
            if frameCount > maxFrames:
                break
            if frame is None:
                break
            Intensity = 0.299 * frame[:,:,0] + 0.587 * frame[:,:,1] + 0.114 * frame[:,:,2]
            # Get size of frame
            height, width = Intensity.shape
            #Green = frame[:,:,1]
            #Blue = frame[:,:,0]
            #newframe = np.zeros((Red.shape[0],Red.shape[1],3,maxFrames+1), np.uint8)
            #newframe[:,:,0,frameCount] = Blue
            #newframe[:,:,1,frameCount] = Green
            #newframe[:,:,2,frameCount] = Red
            if frameCount == 0:
                # Define the codec and create VideoWriter object
                fourcc = cv2.VideoWriter_fourcc(*"MJPG")
                out = cv2.VideoWriter(currFile[:-4] + 'newvid' + '.avi',fourcc, 30.0, (width,height))
            o = out.write(frame.astype('uint8'))
            print(o)
            frameCount = frameCount + 1

        cap.release()
        out.release()
       
convertLowSpeedAVIs()