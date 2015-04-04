#!/usr/bin/env python

import sys
from wsgiref.simple_server import make_server
import cgi
import numpy
from pedees.elliptic2dDriver import Elliptic2dDriver
from pedees.plot import Plot

WIDTH = 300

FORM = b'''
<!DOCTYPE html>

<html>

<head>
<style>
#topPane {
  background-color: lightgray;
  color: black;
  float:left;
  width: %(width)spx;
}
#canvasPane {
  background-color: white;
  color: black;
  width: %(width)spx;
}
#messagePane {
  background-color: black;
  color: red;
  width: %(width)spx;
  text-align: center;
}
</style>
</head>

<title>Solve 2D elliptic equation</title>
</head>
<body>
<form method="post">

<div id="topPane">
<label> Solve - &nabla; 
   <input type="text" placeholder="Enter function(x, y)" name="fFunc" size=10 required>
&nabla; &psi; + 
   <input type="text" placeholder="Enter function(x, y)" name="gFunc" size=10 required>
&psi; =  
   <input type="text" placeholder="Enter function(x, y)" name="sFunc" size=10 required>

<p>

</label> subject to &part; &psi; / &part; n + 
   <input type="text" placeholder="Enter function(t)" name="bFunc" size=13 required>
&psi; =  
   <input type="text" placeholder="Enter function(t)" name="cFunc" size=17 required>
</label> 

</p>

<p>
on x =  <input type="text" placeholder="Enter function(t)" name="xFunc" value="cos(2*pi*t)" size=19 required>
   y = <input type="text" placeholder="Enter function(t)" name="yFunc" value="sin(2*pi*t)" size=19 required>
</label>
<input type="submit" value="Go">
</p></div>
</form>

<p>
<!- the canvas goes here -->
%(canvas)s
</p>

<div id="messagePane">
%(message)s 
</div>
</body>
</html>
'''

def validateXYFunction(str):
  from math import sqrt, log, exp, \
       sin, cos, tan, \
       asin, acos, atan, atan2, \
       pi, sinh, cosh, tanh, asinh, acosh, atanh, e, erf, \
       pow, log10, floor, gamma, factorial, floor, fabs
  x, y = 0., 0.
  msg = ''
  try:
    eval(str)
  except:
    msg = 'Invalid expression %s' % str
  return msg

def application(environ, start_response):

  # no canvas initially
  params = {'width': WIDTH, 'canvas': '', 'message': 'Fill in function fields'}

  if environ['REQUEST_METHOD'] == 'POST':
    post_env = environ.copy()
    post_env['QUERY_STRING'] = ''
    post = cgi.FieldStorage(fp = environ['wsgi.input'], 
                        environ = post_env, 
                        keep_blank_values = True)

    # guard against injection
    fFunc = cgi.escape( post.getfirst("fFunc") )
    gFunc = cgi.escape( post.getfirst("gFunc") )
    sFunc = cgi.escape( post.getfirst("sFunc") )
    xFunc = cgi.escape( post.getfirst("xFunc") )
    yFunc = cgi.escape( post.getfirst("yFunc") )

    # validate
    inputIsValid = True
    for k, func in ('f', fFunc), ('g', gFunc), ('s', sFunc):
      msg = validateXYFunction(func)
      if msg:
        params['message'] = k + ': ' + msg
        inputIsValid = False

    if inputIsValid:
      elliptic = Elliptic2dDriver(fFunc=fFunc, gFunc=gFunc, sFunc=sFunc, \
                                xFunc=xFunc, yFunc=yFunc)

      elliptic.triangulate(numCells = 100)

      bFunc = post.getfirst("bFunc", '1.0')
      cFunc = post.getfirst("cFunc", '1.0')
      elliptic.applyBoundaryConditions(bFunc=bFunc, cFunc=cFunc)

      diag = elliptic.solve()
      params['message'] = 'CG error: ' + diag['error in the norm']

      pl = Plot(elliptic.getTriangulation(), width=WIDTH, height=WIDTH)
      params['canvas'] = '<div id="canvasPane">\n' + pl.jsShow(elliptic.getSolution()) + '\n</div>'

  html = (FORM % params)

  status = '200 OK'
  response_headers = [('Content-Type', 'text/html'),
                      ('Content-Length', str(len(html)))]
  start_response(status, response_headers)
  return [html]

##########################################################################
if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser(description='Run elliptic 2d server.')
  parser.add_argument('--port', dest='port', type=int, 
    default=9000, help='Port number')
  args = parser.parse_args(sys.argv[1:])

  httpd = make_server('localhost', args.port, application)
  httpd.serve_forever()

