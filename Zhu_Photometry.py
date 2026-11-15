import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Zhu_odlib as od
image_data = pd.read_csv('IMAGE_DATA.csv')


exp_time = 90 #seconds
x_values = []
y_values = []
x_day1 = []
y_day1 = []
x_day2 = []
y_day2 = []
x_day3 = []
y_day3 = []
x_day4 = []
y_day4 = []

#day 1
for i in range(2, 7):

    
    #read values from csv
   
    ra1_measured = float(image_data.loc[i, 'ra'])
    dec1_measured = float(image_data.loc[i, 'dec '])
    ra1_ref = float(image_data.loc[i, 'ra (ref)']) / 15
    dec1_ref = float(image_data.loc[i, 'dec (ref)'])
    source_sky = float(image_data.loc[i, 'source-sky'])
    vmag = float(image_data.loc[i, 'vmag (ref)'])
    bmag = float(image_data.loc[i, 'b_mag'])

    #find differences
    diff1_ra = abs(ra1_measured - ra1_ref)
    diff1_dec = abs(dec1_measured - dec1_ref)

    #find RMS error
    RMS1_error_RA = np.sqrt(np.sum(diff1_ra**2)/5)
    RMS1_error_DEC = np.sqrt(np.sum(diff1_dec**2)/5)

    #fit line
    # y1 = vmag - np.log10(source_sky/exp_time)
    # x1 = vmag - bmag
    y1 = vmag
    x1 = np.log10(source_sky/exp_time)
    
    x_values.append(x1)
    y_values.append(y1)
    x_day1.append(x1)
    y_day1.append(y1)

    print("differences1:", diff1_ra, diff1_dec)
    print("RMS_RA error 1:", RMS1_error_RA)
    print("RMS_DEC error 1:", RMS1_error_DEC)

#day 2
for i in range(9, 14):
    #get array
    
    
    ra2_measured = 15*float(image_data.loc[i, 'ra'])
    dec2_measured = float(image_data.loc[i, 'dec '])
    ra2_ref = float(image_data.loc[i, 'ra (ref)'])
    dec2_ref = float(image_data.loc[i, 'dec (ref)'])
    source_sky = float(image_data.loc[i, 'source-sky'])
    vmag = float(image_data.loc[i, 'vmag (ref)'])
    bmag = float(image_data.loc[i, 'b_mag'])

    #find differences
    diff2_ra = abs(ra2_measured - ra2_ref)
    diff2_dec = abs(dec2_measured - dec2_ref)

    #find RMS error
    RMS2_error_RA = np.sqrt(np.sum(diff2_ra**2)/5)
    RMS2_error_DEC = np.sqrt(np.sum(diff2_dec**2)/5)
    #fit line
    # y2 = vmag - np.log10(source_sky/exp_time)
    # x2 = vmag - bmag
    y2 = vmag
    x2 = np.log10(source_sky/exp_time)
    x_values.append(x2)
    y_values.append(y2)
    x_day2.append(x2)
    y_day2.append(y2)
    
    print("differences2:", diff2_ra, diff2_dec)
    print("RMS error 2:", RMS2_error_RA)
    print("RMS_DEC error 2:", RMS2_error_DEC)
    
#day 3
for i in range(16, 21):
    
    
    ra3_measured = 15*float(image_data.loc[i, 'ra'])
    dec3_measured = float(image_data.loc[i, 'dec '])
    ra3_ref = float(image_data.loc[i, 'ra (ref)'])
    dec3_ref = float(image_data.loc[i, 'dec (ref)'])
    source_sky = float(image_data.loc[i, 'source-sky'])
    vmag = float(image_data.loc[i, 'vmag (ref)'])
    bmag = float(image_data.loc[i, 'b_mag'])

    #find differences
    diff3_ra = abs(ra3_measured - ra3_ref)
    diff3_dec = abs(dec3_measured - dec3_ref)

    #find RMS error
    RMS3_error_RA = np.sqrt(np.sum(diff3_ra**2)/5)
    RMS3_error_DEC = np.sqrt(np.sum(diff3_dec**2)/5)
    #fit line
    # y3 = vmag - np.log10(source_sky/exp_time)
    # x3 = vmag - bmag
    y3 = vmag
    x3 = np.log10(source_sky/exp_time)
    x_values.append(x3)
    y_values.append(y3)
    x_day3.append(x3)
    y_day3.append(y3)

    print("differences3:", diff3_ra, diff3_dec)
    print("RMS_RA error 3:", RMS3_error_RA)
    print("RMS_DEC error 3:", RMS3_error_DEC)

#day 4
for i in range(23, 28):
    
    ra4_measured = 15*float(image_data.loc[i, 'ra'])
    dec4_measured = float(image_data.loc[i, 'dec '])
    ra4_ref = float(image_data.loc[i, 'ra (ref)'])
    dec4_ref = float(image_data.loc[i, 'dec (ref)'])
    source_sky = float(image_data.loc[i, 'source-sky'])
    vmag = float(image_data.loc[i, 'vmag (ref)'])
    bmag = float(image_data.loc[i, 'b_mag'])

    #find differences
    diff4_ra = abs(ra4_measured - ra4_ref)
    diff4_dec = abs(dec4_measured - dec4_ref)

    #find RMS error
    RMS4_error_RA = np.sqrt(np.sum(diff4_ra**2)/5)

    #fit line
    # y4 = vmag - np.log10(source_sky/exp_time)
    # x4 = vmag - bmag
    y4 = vmag
    x4 = np.log10(source_sky/exp_time)
    x_values.append(x4)
    y_values.append(y4)
    x_day4.append(x4)
    y_day4.append(y4)

    print("differences4:", diff4_ra, diff4_dec)
    print("RMS error 4:", RMS4_error_RA)


#Make the Plot

# get arrays
x1 = np.array(x_day1)
y1 = np.array(y_day1)

x2 = np.array(x_day2)
y2 = np.array(y_day2)

x3 = np.array(x_day3)
y3 = np.array(y_day3)

x4 = np.array(x_day4)
y4 = np.array(y_day4)

x_total = np.array(x_values)
y_total = np.array(y_values)


#Make plots
a1, b1 = np.polyfit(x1, y1, 1)
plt.scatter(x1, y1, color = 'blue', label = 'Estimating Asteroid Magnitude Night 1')
plt.plot(x1, a1*x1+b1)
plt.show()


a2, b2 = np.polyfit(x2, y2, 1)
plt.scatter(x2, y2, color = 'blue', label = 'Estimating Asteroid Magnitude Night 2')
plt.plot(x2, a2*x2+b2)
plt.show()

a3, b3 = np.polyfit(x3, y3, 1)
plt.scatter(x3, y3, color = 'blue', label = 'Estimating Asteroid Magnitude Night 3')
plt.plot(x3, a3*x3+b3)
plt.show()

a4, b4 = np.polyfit(x4, y4, 1)
plt.scatter(x4, y4, color = 'blue', label = 'Estimating Asteroid Magnitude Night 4')
plt.plot(x4, a4*x4+b4)
plt.show()

a, b = np.polyfit(x_total, y_total, 1)
plt.scatter(x_total, y_total, color = 'blue', label = 'Estimating Asteroid Magnitude')
plt.plot(x_total, a*x_total+b)

#get V mag
V1 = b1 + a1*np.log10(source_sky/exp_time)
V2 = b2 + a2*np.log10(source_sky/exp_time) 
V3 = b3 + a3*np.log10(source_sky/exp_time) 
V4 = b4 + a4*np.log10(source_sky/exp_time) 
V = b + a*np.log10(source_sky/exp_time) 

#Actual V mag
V1_actual = 14.831
V2_actual = 14.921
V3_actual = 15.027
V4_actual = 15.108

#get percent errors
V1_error = 100*abs(V1 - V1_actual)/V1_actual
V2_error = 100*abs(V2 - V2_actual)/V2_actual
V3_error = 100*abs(V3 - V3_actual)/V3_actual
V4_error = 100*abs(V4 - V4_actual)/V4_actual

print("a, b values:")
print("Night 1:", a1, b1)
print("Night 2:", a2, b2)
print("Night 3:", a3, b3)
print("Night 4:", a4, b4)
print("All Nights:", a, b)

print("Apparent Magnitude:")
print("Night 1:", V1)
print("Night 2:", V2)
print("Night 3:", V3)
print("Night 4:", V4)
print("All Nights:", V)

print("Apparent Magnitude % error:")
print("Night 1", V1_error)
print("Night 2", V2_error)
print("Night 3", V3_error)
print("Night 4", V4_error)
