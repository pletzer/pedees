#!/usr/bin/env python

from wsgiref.simple_server import make_server
import cgi
from math import cos, sin, tan, atan, atan2, pi, log, exp, e
import numpy
import random


  def jsDraw(self, triangulation, solution):
    width = 500
    height = int(width * (self.ymax - self.ymin)/(self.xmax - self.xmin))
    res = '''
<canvas id="myCanvas" width="%d" height="%d" style="border:1px solid #000000;">
  <script>
    var c = document.getElementById("myCanvas");
    var ctx = c.getContext("2d");
    <!- fill the triangles -->
''' % (width, height)
  for i in range(len(solution)):
    
    
</canvas>
    '''


def drawBoundary(xExpression, yExpression, numSegments):
  """
  @param xExpression x-expression of 0 <= t <= 1
  @param yExpression y-expression of 0 <= t <= 1
  @param numSegments number of segments
  """
  dt = 1.0 / float(numSegments)
  ts = [i*dt for i in range(numSegments)]
  xs = [eval(xExpression) for t in ts]
  ys = [eval(yExpression) for t in ts]
  xmin = min(xs); ymin = min(ys)
  xmax = max(xs); ymax = max(ys)
  # 1:1 aspect ratio
  height = width * (ymax - ymin) / float(xmax - xmin)
  xPix = int( width*(xs[0] - xmin)/(xmax - xmin) )
  yPix = int( height*(ymax - ys[0])/(ymax - ymin) )
  javaScript = """
<div id="bottomPane">
  var xmin = %f;
  var xmax = %f;
  var ymin = %f;
  var ymax = %f;
  var height = $screen.width (ymax - ymin) / float(xmax - xmin);
  <canvas id="myCanvas" width="$screen.width" height="%d"
  style="border:1px solid #000000;">
    <script>
      var c = document.getElementById("myCanvas");
      var ctx = c.getContext("2d");
      ctx.beginPath();
      ctx.moveTo(%d,%d);
""" % (xmin, xmax, ymin, ymax, xPix, yPix)
  for i in range(1, numSegments + 1):
    xPix = int( width*(xs[i % numSegments] - xmin)/(xmax - xmin) )
    yPix = int( height*(ymax - ys[i % numSegments])/(ymax - ymin) )
    javaScript += """
       ctx.lineTo(%d,%d);
""" % (xPix, yPix)
  javaScript += """
    ctx.stroke();
    </script>
  </canvas>
</div>
"""
  return javaScript

FORM = b'''
<!DOCTYPE html>

<html>

<head>
<style>
#topPane {
	background-color: orange;
	color: black;
	padding:5px;
	float:left;
}
#bottomPane {
	background-color: white;
	color: black;
	float:left;
	padding:5px;
}
</style>
</head>

<title>Solve 2D elliptic equation</title>
</head>
<body>
<form method="post">

<div id="topPane">
<label> Solve - &nabla; 
   <input type="text" placeholder="Enter function(x, y)" name="fFunc" value="1" size=10>
&nabla; &psi; + 
   <input type="text" placeholder="Enter function(x, y)" name="gFunc" value="0" size=10>
&psi; =  
   <input type="text" placeholder="Enter function(x, y)" name="sFunc" value="0" size=10>

</label> subject to <input type="text" placeholder="Enter function(t)" name="aFunc" value="0" size=2> 
&part; &psi; / &part; n + 
   <input type="text" placeholder="Enter function(t)" name="bFunc" value="1" size=2>
&psi; =  
   <input type="text" placeholder="Enter function(t)" name="cFunc" value="0" size=2>
</label> 
on x =  <input type="text" placeholder="Enter function(t)" name="xFunc" value="cos(2*pi*t)" size=10>
   y = <input type="text" placeholder="Enter function(t)" name="yFunc" value="sin(2*pi*t)" size=10>
</label>
<input type="submit" value="Go">
</div>
</form>

<!- the canvas goes here -->
%s
</body>
</html>
'''

def application(environ, start_response):

	html = (FORM % '')

	if environ['REQUEST_METHOD'] == 'POST':
		post_env = environ.copy()
		post_env['QUERY_STRING'] = ''
		post = cgi.FieldStorage(fp = environ['wsgi.input'], 
				                environ = post_env, 
				                keep_blank_values = True)

    elliptic = Elliptic2D(f=post['fFunc'], 
                          g=post['gFunc'],
                          s=post['sFunc'],
                          xBound=post['xFunc'],
                          yBound=post['yFunc'],
                          bcs=post['BCs'])

    js = draw(elliptic, width=500)

		boundaryPointsJS = drawBoundary(width=500, 
      xExpression=post['xFunc'].value, 
      yExpression=post['yFunc'].value, 
      numSegments=40)
		html = (FORM % boundaryPointsJS)

	status = '200 OK'
	response_headers = [('Content-Type', 'text/html'),
	                    ('Content-Length', str(len(html)))]
	start_response(status, response_headers)
	return [html]

##########################################################################
if __name__ == '__main__':

	import sys
	import argparse

	parser = argparse.ArgumentParser(description='Run elliptic 2d server.')
	parser.add_argument('--port', dest='port', type=int, 
    default=9000, help='Port number')
	args = parser.parse_args(sys.argv[1:])

	httpd = make_server('localhost', args.port, application)
	httpd.serve_forever()

