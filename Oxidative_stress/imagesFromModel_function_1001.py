# -*- coding: utf-8 -*-
"""
@author: some4879
"""
from PIL import Image, ImageDraw, ImageFilter
import random
import matplotlib.pyplot as plt
import tifffile
import os

list2=[0.] 
x = 0;y = 0+40;dx = 10 ;dy =11+40  
class main():
    
    def make_ellipse(self, image, x, y, dx, dy, next_1, color):
        next_1=next_1-10
        dr = ImageDraw.Draw(image)
        ell = dr.ellipse((x+1+2, y, dx+1+2, dy-1), fill=color)
        ell2 = dr.ellipse((x+next_1+2 , y, dx+next_1+2, dy-1), color)
        img1 = Image.new("RGB", (next_1, dy - y), color=color)## makes the rectangle
        image.paste(img1, ((  int((dx)), int(y))))## where rectangle paste with y telling the rightmost point
        return ell, ell2, img1

    def draw_ellipses(self,length, intt):
        img = Image.new("RGBA", (160, 160), "rgba(0,0,0,0)") 
        ell1 = main.make_ellipse(self,img, x , y , dx , dy , int(length-1 ), "rgb(%d,%d,%d)"%(intt,intt,intt)) 
        return img
 
    def image_edit2(self,img):
        img = img.rotate(0)
        imgSmall = img.resize((600, 600), resample=Image.BILINEAR)
        result = imgSmall.resize(img.size, Image.NEAREST)
        return result

    def tps2(self,ll,timedur , trenches, cellnn  ,gett__,cellnn_ ,  alllgrowths_,grxaall,h2o2_,oxyRall,katgall,ahpcall,numberoftrench ) :
          timedurarr=[]
          lenarr=[]

          isExist = os.path.exists('Output_images')
          if not isExist:
               os.makedirs('Output_images')

          for l in range(0, timedur-2):
             img1 = Image.new("RGB", (50*5*2 , 60*(numberoftrench)), "rgba(18,18,18,100)")
             for trenchNum in range(len(trenches)):
                 tt = trenches[trenchNum]
                 for cc in range(len(cellnn_[l][trenchNum])):   
                    length1 = gett__[l][trenchNum][cc]; totaa = gett__[l][trenchNum][0:cc]; tota = sum ( totaa )
                    timedurarr+=[l]; lenarr+=[int(tota)   ] 
                    img = main.draw_ellipses( self,int(length1),  40+int( grxaall[l] [trenchNum][cc]   *50)  ) 
                    result =img
                    img1.paste(result, (int(tota),trenchNum*50), result)   ##where new cell pasted
             img1 = main.image_edit2(self,img1)
             img2 = img1.rotate(-90, expand=True)    
             plt.imshow(img2)
             plt.title('%s'%l)
             resized_img = img2.resize((60*(numberoftrench)*1, 50*5*2 )); tifffile.imsave("Output_images/movie-xy%d.tiff" % ( l), resized_img)
          return  cellnn  ,   length1 ,timedurarr ,lenarr#'''
