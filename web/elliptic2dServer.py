#!/usr/bin/env python

from wsgiref.simple_server import make_server
import cgi
import numpy
from pedees.ellipt2dDriver import Ellipt2dDriver
from pedees import plot


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

    elliptic = Elliptic2DDriver(fFunc=post['fFunc'], 
                          gFunc=post['gFunc'],
                          sFunc=post['sFunc'],
                          xFunc=post['xFunc'],
                          yFunc=post['yFunc'])
    elliptic.triangulate(numCells = 100)

    elliptic.applyBoundaryConditions(bFunc=post['bFunc'], cFunc=post['cFunc'])

    elliptic.solve()

    js = plot.jsShow(elliptic, width=500, height=500)

		html = (FORM % js)

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

