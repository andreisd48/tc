import numpy as np
import collections
import matplotlib.pyplot as plt
import h5py

number_channels = 64
thetaDivider = 50
hdf5 = 'DeepInsighthdf5.h5'
target_frequency = 850
tf_list = np.zeros(thetaDivider, number_channels)
cf_list = np.zeros(thetaDivider, number_channels)
tf_time_list = np.array(thetaDivider,2)
cf_time_list = np.zeros(thetaDivider,2)
wavelets = np.array() #should be a dataset of hdf5 when implemented
fourier_frequencies = np.array() #should be a dataset of hdf5 when implemented
head_direction = np.array() #should be a dataset of hdf5 when implemented
tf_angle_counter = np.zeros(thetaDivider)
cf_angle_counter = np.zeros(thetaDivider)
final_list = collections.OrderedDict

if (target_frequency < 5000):
    control_frequency = target_frequency * 2
else:
    control_frequency = target_frequency / 2

#hdf5_file = h5py.File(hdf5, mode='r+')
#fourier_frequencies = hdf5_file["inputs/fourier_frequencies"][()].astype(np.float32)

def getAngles(thetaDivider):
    angles = np.linspace(-np.pi,np.pi,thetaDivider)
    return angles

def freq_index(fourier_frequencies, target_frequency, control_frequency):
    for i in fourier_frequencies:
        if (i == target_frequency):
            tf = i
        if (i == control_frequency):
            cf = i
    return tf, cf

def tuningCurve(hdf5_file, channel_number, target_frequency, control_frequency):
    #separate the hdf5_file variable into its parts: inputs/wavelets and outputs/head_direction
    tf_index, cf_index = freq_index(fourier_frequencies, target_frequency, control_frequency) #get the index of the target and control index of the 25 frequencies
    i = 0
    angle_index = 0

    # go through every timestep in the hdf5 file, and if the head direction at that time matches one of the theta values we're interested in, add that amplitude to the
    #array for that given frequency or the control (first if statement). Add one to the counter for the number of frequencies found at that angle, iterate to next timestep

    while (i < hdf5_file.size[0]):
        if (head_direction[i] in getAngles(thetaDivider)): #get index so you know which angle it matched with
            angle_index = np.where(head_direction[i] == getAngles(thetaDivider))
            if (head_direction[i] == target_frequency): #make sure this is inputs/wavelets when actually running
                tf_list[angle_index] += hdf5_file[i,tf_index,:]
                tf_angle_counter[angle_index] += 1
                i+= 1
                if (tf_time_list[angle_index,0] == 0): #if the time stamp currently is 0 then make the 1st index the current time, if its not make second index current time
                    tf_time_list[angle_index, 0] == hdf5_file[i,0,0]
                else:
                    tf_time_list[angle_index, 1] == hdf5_file[i,0,0]
            if (head_direction[i] == control_frequency):
                cf_list[angle_index] += hdf5_file[i,cf_index,:]
                cf_angle_counter[angle_index] += 1
                i+= 1
                if (cf_time_list[angle_index,0] == 0):
                    cf_time_list[angle_index, 0] == hdf5_file[i,0,0]
                else:
                    cf_time_list[angle_index, 0] == hdf5_file[i,0,0]
        else:
            i += 1

    #call function to calculate the tuning curve per angle. Already have the amplitudes added up, just need to divde by timstep difference and then # of frequencies found
    #at that angle

    for i in range(thetaDivider):
        tf_list[i] = (tf_list[i])/(tf_angle_counter[i] * (tf_time_list[i,1] - tf_time_list[i,0]))
        cf_list[i] = (cf_list[i])/(cf_angle_counter[i] * (cf_time_list[i,1] - cf_time_list[i,0]))

    #using matplotlib, plot 128 tuning curves, one for each of the 64 channel, with theta_value # of datapoints, one for each angle for control and target frequency

    for i in range(channel_number):
        tf_x = getAngles(thetaDivider)
        cf_x = getAngles(thetaDivider)
        tf_y = tf_list[:,i]
        cf_y = cf_list[:,i]

        plt.plot(tf_x,tf_y)
        plt.xlabel("theta")
        plt.ylabel("average amplitude")
        plt.title("tf: average amplitude per theta for channel ", i)
        plt.show()

        plt.plot(cf_x, cf_y)
        plt.xlabel("theta")
        plt.ylabel("average amplitude")
        plt.title("cf: average amplitude per theta for channel ", i)
        plt.show()


