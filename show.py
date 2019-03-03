import cv2
import numpy as np
from gifwriter import GifWriter
select='0'
size=500
gif_writer = GifWriter('temp/%03d.png', 'output/process_korea'+select+'.gif')
green = (0, 255, 0)
city_color=[(150,150,255),(0,0,255),(0,0,0)]
f_process=open('data/process_korea'+select+'.dat', 'r')
for iter in range(1,int(f_process.readline().split()[1])+1):
    img = 255 * np.ones((size+50, size, 3), dtype=np.uint8)
    f_loc= open("data/loc.dat", 'r')
    f_ts=open("data/ts.dat",'r')
    Q=[float(s) for s in f_process.readline().split()]
    for i in range(int(f_ts.readline())):
        poly=[]
        for j in range(int(f_ts.readline().split()[1])):
            poly.append([int(float(s)*size) for s in f_ts.readline().split()])
        poly=np.array(poly,np.int32)
        cv2.fillPoly(img, [poly], (255,255*(1-Q[i]),255*(1-Q[i])))
        cv2.polylines(img, [poly], True, green, 1)
    for i in range(int(f_loc.readline())):
        data=f_loc.readline().split()
        cv2.circle(img,tuple([(int)(size*float(s)) for s in data[1:3]]), 3,city_color[int(data[0])-1],-1)
    f_loc.close()
    f_ts.close()
    cv2.putText(img, 'korea '+select+': t='+str(iter), (size>>2, size+30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0,0,0), 2)
    gif_writer.append(img)
f_process.close()
gif_writer.close()


