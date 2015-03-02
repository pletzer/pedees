#!/usr/bin/env python

import Tkinter
import math

class Plot:

  def __init__(self, tri, width, height):

    self.padding = 20
    self.w = width - 2*self.padding
    self.h = height - 2*self.padding

    self.master = Tkinter.Tk()
    self.canvas = Tkinter.Canvas(bg="white", width=width, height=height)
    self.canvas.pack()

    self.triangulation = tri
    points = tri.getPoints()
    self.xmin = min( [p[0] for p in points] )
    self.xmax = max( [p[0] for p in points] )
    self.ymin = min( [p[1] for p in points] )
    self.ymax = max( [p[1] for p in points] )

  def x2Pix(self, x):
    return self.padding + self.w*(x - self.xmin)/(self.xmax - self.xmin)

  def y2Pix(self, y):
    return self.padding + self.h*(self.ymax - y)/(self.ymax - self.ymin)

  def getRGB(self, v, vmin, vmax):
    x = (v - vmin)/(vmax - vmin)
    blue = min((max((4*(0.75-x), 0.)), 1.))
    red  = min((max((4*(x-0.25), 0.)), 1.))
    green= min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
    return "#%02x%02x%02x" % (red*255, green*255, blue*255)

  def jsShow(self, solution, width, height):
    res = '''
<canvas id="myCanvas" width="%d" height="%d" style="border:1px solid #000000;">
  <script>
    var c = document.getElementById("myCanvas");
    var ctx = c.getContext("2d");
    <!- fill the triangles -->
          ''' % (width, height)

    vmin, vmax = min(solution), max(solution)
    n = len(solution)
    points = self.triangulation.getPoints()

    for cell in self.triangulation.getTriangles():
      ia, ib, ic = cell
      va, vb, vc = points[ia], points[ib], points[ic]
      va, vb, vc = points[ia], points[ib], points[ic]
      pxa, pxb, pxc = self.x2Pix(va[0]), self.x2Pix(vb[0]), self.x2Pix(vc[0])
      pya, pyb, pyc = self.y2Pix(va[1]), self.y2Pix(vb[1]), self.y2Pix(vc[1])
      val = (solution[ia] + solution[ib] + solution[ic])/3.0
      r, g, b = self.getRGB(val, vmin, vmax)
      colorHexStr = '%x' % r, '%x' % g, '%x' % b
      res += '''
        ctx.beginPath(); ctx.moveTo(%d, %d); ctx.lineTo(%d, %d); ctx.lineTo(%d, %d); 
        ctx.fillStyle(%s); ctx.fill();
        ''' % (pxa, pya, pxb, pyb, pxc, pyc, colorHexStr)
      res += '''    
  </script> 
</canvas>
      '''
    return res

  def show(self, solution):

    vmin, vmax = min(solution), max(solution)
    n = len(solution)
    points = self.triangulation.getPoints()

    for cell in self.triangulation.getTriangles():
      ia, ib, ic = cell
      va, vb, vc = points[ia], points[ib], points[ic]
      pxa, pxb, pxc = self.x2Pix(va[0]), self.x2Pix(vb[0]), self.x2Pix(vc[0])
      pya, pyb, pyc = self.y2Pix(va[1]), self.y2Pix(vb[1]), self.y2Pix(vc[1])
      val = (solution[ia] + solution[ib] + solution[ic])/3.0
      color = self.getRGB(val, vmin, vmax)
      self.canvas.create_polygon(pxa, pya, pxb, pyb, pxc, pyc, fill=color)

    for edge in self.triangulation.getEdges():
      i1, i2 = edge
      v1, v2 = points[i1], points[i2]
      px1, px2 = self.x2Pix(v1[0]), self.x2Pix(v2[0])
      py1, py2 = self.y2Pix(v1[1]), self.y2Pix(v2[1])
      self.canvas.create_line(px1, py1, px2, py2, fill='black')

    for i in range(len(points)):
      v = points[i]
      px, py = self.x2Pix(v[0]) - 5, self.y2Pix(v[1]) - 5
      self.canvas.create_text(px, py, text=str(i))
      px += 10
      py += 10
      self.canvas.create_text(px, py, text='%5.2f'%solution[i])

    self.master.mainloop()
