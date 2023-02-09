# for all .AVI files in directory, convert to smaller .avi files with different codec

import numpy as np
import cv2
import os
import sys
import scipy.io
#import matplotlib.pyplot as plt

def convertLowSpeedAVIs():

    # n_frames_before_event
    n_frames_before_event = 4
    # n_frames_after_event
    n_frames_after_event = 196

    i=0
    useTheseFiles = []
    useTheseFolders = []
    dirname = sys.argv[1]
    for file in sorted(os.listdir(dirname)):     
        if file.endswith(".AVI") | file.endswith(".avi"):
            print(os.path.join(dirname, file))
            useTheseFiles.append(os.path.join(dirname, file))
            if file.endswith(".avi"):
                thisisshortvid = True
                #useTheseFolders.append(os.path.join(dirname, file.split('_')[0]) + '*' + '_processed_data')
                useTheseFolders.append(file.split('_')[0] + '*' + '_processed_data')
            else:
                thisisshortvid = False
                # get associated processed_data folder
                #useTheseFolders.append(os.path.join(dirname, file)[:-4] + '_processed_data')
                useTheseFolders.append(file[:-4] + '_processed_data')
            i=i+1

    for x in range(i):
        currFile = useTheseFiles[x]
        print("Now processing: " + currFile)
        
        cap = cv2.VideoCapture(currFile)

        maxFrames = 50000
        frameCount = 0
        whichsmallvid = 0

        if thisisshortvid == True:
            # If variable tbt does not yet exist
            if 'tbt' not in locals():
                if '*' in useTheseFolders[x]:
                    for folder in sorted(os.listdir(dirname)):
                        if 'processed_data' in folder:
                            tbt = scipy.io.loadmat(os.path.join(dirname, folder) + '/tbt.mat')
                            break
                else:
                    tbt = scipy.io.loadmat(useTheseFolders[x] + '/tbt.mat')
                # index field in matlab struct
                tbt = tbt['tbt'][0,0]
                #temp = tbt['reachBatch_success_reachStarts']
                #temp = temp.flatten()
                # plot temp
                #plt.plot(temp)
                #plt.show()
                # get behavior events from tbt
                eventsInMovieFrames = getBehaviorEvents(tbt, n_frames_before_event, n_frames_after_event)
                for key in eventsInMovieFrames:
                    print(eventsInMovieFrames[key])
                    print(2481 in eventsInMovieFrames[key])
            else:
                # Adjust eventsInMovieFrames to account for short vid
                for key in eventsInMovieFrames:
                    eventsInMovieFrames[key] = [x - shortvidframecount for x in eventsInMovieFrames[key]]
        else:
            # read in tbt.mat from processed_data folder
            tbt = scipy.io.loadmat(useTheseFolders[x] + '/tbt.mat')
            # index field in matlab struct
            tbt = tbt['tbt'][0,0]
            #temp = tbt['reachBatch_success_reachStarts']
            #temp = temp.flatten()
            # plot temp
            #plt.plot(temp)
            #plt.show()
            # get behavior events from tbt
            eventsInMovieFrames = getBehaviorEvents(tbt, n_frames_before_event, n_frames_after_event)
            for key in eventsInMovieFrames:
                print(eventsInMovieFrames[key])
                print(2481 in eventsInMovieFrames[key])

        while(cap.isOpened()):
            ret, frame = cap.read()
            # Start frame count for short vids
            shortvidframecount = 0
            if frameCount%500 == 0:
                print(frameCount)
            if frameCount > maxFrames:
                break
            if frame is None:
                break
            # Get size of frame
            height, width = frame[:,:,0].shape
            if frameCount%5000 == 0:
                if whichsmallvid != 0:
                    out.release()
                    f.close()
                # Define the codec and create VideoWriter object
                fourcc = cv2.VideoWriter_fourcc(*"MJPG")
                # format whichsmallvid to 4 digit string
                out = cv2.VideoWriter(currFile[:-4] + 'newvid' + str(whichsmallvid).zfill(4) + '.avi',fourcc, 30.0, (width,height))
                # open a .csv file to write frame numbers
                f = open(currFile[:-4] + 'newvid' + str(whichsmallvid).zfill(4) + '.csv', 'w')
                smallFileFrameCount = 1
                # write header to .csv file
                f.write('Frame, success, drop, miss, missing_pellet\n')
                whichsmallvid = whichsmallvid + 1
            out.write(frame.astype('uint8'))
            # write frame number to first column of .csv file
            f.write(str(smallFileFrameCount) + ',')
            # write eventsInMovieFrames as next columns of .csv file
            for key in eventsInMovieFrames:
                # if frameCount is in eventsInMovieFrames[key], write 1, else write 0
                if frameCount in eventsInMovieFrames[key]:
                    f.write('1,')
                else:
                    f.write('0,')
            # write newline to .csv file
            f.write('\n')
            frameCount = frameCount + 1
            smallFileFrameCount = smallFileFrameCount + 1
            shortvidframecount = shortvidframecount + 1

        cap.release()
        

def getBehaviorEvents(tbt, n_frames_before_event, n_frames_after_event):
    movieframeinds = tbt['movieframeinds']
    # turn 2D movieframeinds into 1D array
    movieframeinds = movieframeinds.flatten()
    # make a dictionary of behavior events called eventsInMovieFrames
    eventsInMovieFrames = {}
    eventsInMovieFrames['success']=getFramesForBehEvents(movieframeinds, tbt['reachBatch_success_reachStarts'] + tbt['reachBatch_success_reachStarts_pawOnWheel'], n_frames_before_event, n_frames_after_event)
    eventsInMovieFrames['drop']=getFramesForBehEvents(movieframeinds, tbt['reachBatch_drop_reachStarts'] + tbt['reachBatch_drop_reachStarts_pawOnWheel'], n_frames_before_event, n_frames_after_event)
    eventsInMovieFrames['miss']=getFramesForBehEvents(movieframeinds, tbt['reachBatch_miss_reachStarts'] + tbt['reachBatch_miss_reachStarts_pawOnWheel'], n_frames_before_event, n_frames_after_event)
    eventsInMovieFrames['missing_pellet']=getFramesForBehEvents(movieframeinds, tbt['pelletmissingreach_reachStarts'], n_frames_before_event, n_frames_after_event)
    return eventsInMovieFrames

def getFramesForBehEvents(movieframeinds, behaviorEvents, n_frames_before_event, n_frames_after_event):
    behaviorEvents = behaviorEvents.flatten()
    behaviorEvents[behaviorEvents>1] = 1
    # get values of movieframeinds where behaviorEvents are 1
    mvinds = movieframeinds[behaviorEvents==1]
    # get floor of mvinds
    mvinds = np.floor(mvinds)
    # add to mv_inds integers in the range from mvinds-n_frames_before_event to mvinds+n_frames_after_event
    temp = [range(x-n_frames_before_event,x+n_frames_after_event) for x in mvinds.astype(int)]
    temp=np.array(temp).flatten()
    # add temp to mvinds
    mvinds = np.concatenate((mvinds, temp))
    mvinds[mvinds<1] = 1
    mvinds[mvinds>movieframeinds[-1]] = movieframeinds[-1]
    # get unique mvinds
    mvinds = np.unique(mvinds)
    # convert to int
    mvinds = mvinds.astype(int)
    return mvinds




convertLowSpeedAVIs()