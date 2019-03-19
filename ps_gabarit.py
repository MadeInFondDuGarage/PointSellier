#!/usr/bin/env python
'''
Copyright (C) 2005 Aaron Spike, aaron@ekips.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
'''
import inkex, simplestyle, simplepath, math
import cubicsuperpath
import bezmisc
import re
import cmath
import pathmodifier
import simpletransform

class P:
	def __init__(self, x, y=None):
		if not y==None:
			self.x, self.y = float(x), float(y)
		else:
			self.x, self.y = float(x[0]), float(x[1])
	def __add__(self, other): return P(self.x + other.x, self.y + other.y)
	def __sub__(self, other): return P(self.x - other.x, self.y - other.y)
	def __neg__(self): return P(-self.x, -self.y)
	def __mul__(self, other):
		if isinstance(other, P):
			return self.x * other.x + self.y * other.y
		return P(self.x * other, self.y * other)
	__rmul__ = __mul__
	def __div__(self, other): return P(self.x / other, self.y / other)
	def mag(self): return math.hypot(self.x, self.y)
	def unit(self):
		h = self.mag()
		if h: return self / h
		else: return P(0,0)
	def dot(self, other): return self.x * other.x + self.y * other.y
	def rot(self, theta):
		c = math.cos(theta)
		s = math.sin(theta)
		return P(self.x * c - self.y * s,  self.x * s + self.y * c)
	def angle(self): return math.atan2(self.y, self.x)
	def __repr__(self): return '%f,%f' % (self.x, self.y)
	def pr(self): return "%.2f,%.2f" % (self.x, self.y)
	def to_list(self): return [self.x, self.y]
	def ccw(self): return P(-self.y,self.x)
	def l2(self): return self.x*self.x + self.y*self.y

def csp_subpath_ccw(subpath):
	# Remove all zerro length segments
	s = 0
	#subpath = subpath[:]
	if (P(subpath[-1][1])-P(subpath[0][1])).l2() > 1e-10 :
		subpath[-1][2] = subpath[-1][1]
		subpath[0][0] = subpath[0][1]
		subpath += [ [subpath[0][1],subpath[0][1],subpath[0][1]]  ]
	pl = subpath[-1][2]
	for sp1 in subpath:
		for p in sp1 :
			s += (p[0]-pl[0])*(p[1]+pl[1])
			pl = p
	return s<0

def cubic_solver(a,b,c,d):
	if a!=0:
		#	Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
		a,b,c = (b/a, c/a, d/a)
		m = 2*a**3 - 9*a*b + 27*c
		k = a**2 - 3*b
		n = m**2 - 4*k**3
		w1 = -.5 + .5*cmath.sqrt(3)*1j
		w2 = -.5 - .5*cmath.sqrt(3)*1j
		if n>=0 :
			t = m+math.sqrt(n)
			m1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
			t = m-math.sqrt(n)
			n1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
		else :
			m1 = pow(complex((m+cmath.sqrt(n))/2),1./3)
			n1 = pow(complex((m-cmath.sqrt(n))/2),1./3)
		x1 = -1./3 * (a + m1 + n1)
		x2 = -1./3 * (a + w1*m1 + w2*n1)
		x3 = -1./3 * (a + w2*m1 + w1*n1)
		return [x1,x2,x3]
	elif b!=0:
		det = c**2-4*b*d
		if det>0 :
			return [(-c+math.sqrt(det))/(2*b),(-c-math.sqrt(det))/(2*b)]
		elif d == 0 :
			return [-c/(b*b)]
		else :
			return [(-c+cmath.sqrt(det))/(2*b),(-c-cmath.sqrt(det))/(2*b)]
	elif c!=0 :
		return [-d/c]
	else : return []

def csp_true_bounds(csp) :
	# Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
	minx = [float("inf"), 0, 0, 0]
	maxx = [float("-inf"), 0, 0, 0]
	miny = [float("inf"), 0, 0, 0]
	maxy = [float("-inf"), 0, 0, 0]
	for i in range(len(csp)):
		for j in range(1,len(csp[i])):
			ax,ay,bx,by,cx,cy,x0,y0 = bezmisc.bezierparameterize((csp[i][j-1][1],csp[i][j-1][2],csp[i][j][0],csp[i][j][1]))
			roots = cubic_solver(0, 3*ax, 2*bx, cx)	 + [0,1]
			for root in roots :
				if type(root) is complex and abs(root.imag)<1e-10:
					root = root.real
				if type(root) is not complex and 0<=root<=1:
					y = ay*(root**3)+by*(root**2)+cy*root+y0
					x = ax*(root**3)+bx*(root**2)+cx*root+x0
					maxx = max([x,y,i,j,root],maxx)
					minx = min([x,y,i,j,root],minx)

			roots = cubic_solver(0, 3*ay, 2*by, cy)	 + [0,1]
			for root in roots :
				if type(root) is complex and root.imag==0:
					root = root.real
				if type(root) is not complex and 0<=root<=1:
					y = ay*(root**3)+by*(root**2)+cy*root+y0
					x = ax*(root**3)+bx*(root**2)+cx*root+x0
					maxy = max([y,x,i,j,root],maxy)
					miny = min([y,x,i,j,root],miny)
	maxy[0],maxy[1] = maxy[1],maxy[0]
	miny[0],miny[1] = miny[1],miny[0]

	return minx,miny,maxx,maxy

def csp_split(sp1,sp2,t=.5) :
	[x1,y1],[x2,y2],[x3,y3],[x4,y4] = sp1[1], sp1[2], sp2[0], sp2[1]
	x12 = x1+(x2-x1)*t
	y12 = y1+(y2-y1)*t
	x23 = x2+(x3-x2)*t
	y23 = y2+(y3-y2)*t
	x34 = x3+(x4-x3)*t
	y34 = y3+(y4-y3)*t
	x1223 = x12+(x23-x12)*t
	y1223 = y12+(y23-y12)*t
	x2334 = x23+(x34-x23)*t
	y2334 = y23+(y34-y23)*t
	x = x1223+(x2334-x1223)*t
	y = y1223+(y2334-y1223)*t
	return [sp1[0],sp1[1],[x12,y12]], [[x1223,y1223],[x,y],[x2334,y2334]], [[x34,y34],sp2[1],sp2[2]]

class path_modif(pathmodifier.Diffeo):
    def __init__(self):
        pathmodifier.Diffeo.__init__(self)

    def object_to_path(self,node):
        self.objectToPath(node)

class Dots(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-d", "--dotsize",
                        action="store", type="string",
                        dest="dotsize", default="10mm",
                        help="Size of the dots placed at path nodes")
        self.OptionParser.add_option("--stylegab",
                                     action="store", type="string",
                                     dest="stylegab", default="Gabcarrer",
                                     help="Type de Gabarit")
        self.OptionParser.add_option("--tab",
                        action="store", type="string",
                        dest="tab",
                        help="The selected UI-tab when OK was pressed")

    def separateLastAndFirst(self, p):
        # Separate the last and first dot if they are togheter
        lastDot = -1
        if p[lastDot][1] == []: lastDot = -2
        if round(p[lastDot][1][-2]) == round(p[0][1][-2]) and \
                round(p[lastDot][1][-1]) == round(p[0][1][-1]):
                x1 = p[lastDot][1][-2]
                y1 = p[lastDot][1][-1]
                x2 = p[lastDot-1][1][-2]
                y2 = p[lastDot-1][1][-1]
                dx = abs( max(x1,x2) - min(x1,x2) )
                dy = abs( max(y1,y2) - min(y1,y2) )
                dist = math.sqrt( dx**2 + dy**2 )
                x = dx/dist
                y = dy/dist
                if x1 > x2: x *= -1
                if y1 > y2: y *= -1
                p[lastDot][1][-2] += x
                p[lastDot][1][-1] += y

    def addDot(self, node, carrer_cercle):
        self.dotGroup = inkex.etree.SubElement( node.getparent(), inkex.addNS('g','svg') )
        self.dotGroup.set(inkex.addNS('label','inkscape'), 'gabarit: ' + self.options.dotsize +" " + carrer_cercle)
        try:
            t = node.get('transform')
            self.dotGroup.set('transform', t)
        except:
            pass
        style = simplestyle.formatStyle({ 'stroke': '#00ff00', 'fill': 'none' ,'stroke-width': str(self.unittouu('1px'))})
        a = []
        p = simplepath.parsePath(node.get('d'))
        for cmd,params in p:
            if cmd != 'Z' and cmd != 'z':
                if carrer_cercle=="Carrer":
                    dot_att = {
                      'style': style,
                      'width': str( self.unittouu(self.options.dotsize)*2 ),
                      'height': str( self.unittouu(self.options.dotsize)*2 ),
                      'x': str( params[-2]-self.unittouu(self.options.dotsize) ),
                      'y': str( params[-1] -self.unittouu(self.options.dotsize) )
                    }

                    inkex.etree.SubElement(
                      self.dotGroup,
                      inkex.addNS('rect','svg'),
                      dot_att )
                else:
                    dot_att = {
                      'style': style,
                      'r':  str( self.unittouu(self.options.dotsize) ),
                      'cx': str( params[-2] ),
                      'cy': str( params[-1] )
                    }
                    inkex.etree.SubElement(self.dotGroup,inkex.addNS('circle','svg'),dot_att )

	def transform(self,source_point, layer, reverse=False):
		if layer not in self.transform_matrix:
			for i in range(self.layers.index(layer),-1,-1):
				if self.layers[i] in self.orientation_points :
					break
			if self.layers[i] not in self.orientation_points :
				self.error(_("Orientation points for '%s' layer have not been found! Please add orientation points using Orientation tab!") % layer.get(inkex.addNS('label','inkscape')),"no_orientation_points")
			elif self.layers[i] in self.transform_matrix :
				self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
				self.Zcoordinates[layer] = self.Zcoordinates[self.layers[i]]
			else :
				orientation_layer = self.layers[i]
				if len(self.orientation_points[orientation_layer])>1 :
					self.error(_("There are more than one orientation point groups in '%s' layer") % orientation_layer.get(inkex.addNS('label','inkscape')),"more_than_one_orientation_point_groups")
				points = self.orientation_points[orientation_layer][0]
				if len(points)==2:
					points += [ [ [(points[1][0][1]-points[0][0][1])+points[0][0][0], -(points[1][0][0]-points[0][0][0])+points[0][0][1]], [-(points[1][1][1]-points[0][1][1])+points[0][1][0], points[1][1][0]-points[0][1][0]+points[0][1][1]] ] ]
				if len(points)==3:
					print_("Layer '%s' Orientation points: " % orientation_layer.get(inkex.addNS('label','inkscape')))
					for point in points:
						print_(point)
					#	Zcoordinates definition taken from Orientatnion point 1 and 2
					self.Zcoordinates[layer] = [max(points[0][1][2],points[1][1][2]), min(points[0][1][2],points[1][1][2])]
					matrix = numpy.array([
								[points[0][0][0], points[0][0][1], 1, 0, 0, 0, 0, 0, 0],
								[0, 0, 0, points[0][0][0], points[0][0][1], 1, 0, 0, 0],
								[0, 0, 0, 0, 0, 0, points[0][0][0], points[0][0][1], 1],
								[points[1][0][0], points[1][0][1], 1, 0, 0, 0, 0, 0, 0],
								[0, 0, 0, points[1][0][0], points[1][0][1], 1, 0, 0, 0],
								[0, 0, 0, 0, 0, 0, points[1][0][0], points[1][0][1], 1],
								[points[2][0][0], points[2][0][1], 1, 0, 0, 0, 0, 0, 0],
								[0, 0, 0, points[2][0][0], points[2][0][1], 1, 0, 0, 0],
								[0, 0, 0, 0, 0, 0, points[2][0][0], points[2][0][1], 1]
							])

					if numpy.linalg.det(matrix)!=0 :
						m = numpy.linalg.solve(matrix,
							numpy.array(
								[[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]], [points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
										)
							).tolist()
						self.transform_matrix[layer] = [[m[j*3+i][0] for i in range(3)] for j in range(3)]

					else :
						self.error(_("Orientation points are wrong! (if there are two orientation points they should not be the same. If there are three orientation points they should not be in a straight line.)"),"wrong_orientation_points")
				else :
					self.error(_("Orientation points are wrong! (if there are two orientation points they should not be the same. If there are three orientation points they should not be in a straight line.)"),"wrong_orientation_points")

			self.transform_matrix_reverse[layer] = numpy.linalg.inv(self.transform_matrix[layer]).tolist()

			###self.Zauto_scale[layer]  = math.sqrt( (self.transform_matrix[layer][0][0]**2 + self.transform_matrix[layer][1][1]**2)/2 )
			### Zautoscale is absolete
			self.Zauto_scale[layer] = 1
		x,y = source_point[0], source_point[1]
		if not reverse :
			t = self.transform_matrix[layer]
		else :
			t = self.transform_matrix_reverse[layer]
		return [t[0][0]*x+t[0][1]*y+t[0][2], t[1][0]*x+t[1][1]*y+t[1][2]]

    def effect(self):
        if len(self.options.ids)<=0:
            inkex.errormsg("This extension requires at least one selected path.")
            return
        if self.options.stylegab == "Carrer" or self.options.stylegab =="Cercle":
            selection = self.selected
            if (selection):
                for id, node in selection.iteritems():
                    if node.tag == inkex.addNS('path','svg'):
                        self.addDot(node, self.options.stylegab)
                    else:
                        inkex.errormsg("This extension need a path, not groups.")
        else:
            selection = self.selected
            if (selection):
                for id, path in selection.iteritems():#test pour verifier queje suis sur un chemin et pas un groups,ajout d'un message d'alarme
                    if path.tag == inkex.addNS('path','svg'):
                        if self.options.stylegab == "Interieur":
                            in_out=-1
                        else:
                            in_out=1
                        d = path.get('d')
                        csp = cubicsuperpath.parsePath(d)

                        if path.get(inkex.addNS('type','sodipodi'))!="inkscape:offset":
                            min_x,min_y,min_i,min_j,min_t = csp_true_bounds(csp)[1]

                            if min_y!=float("-inf") :
                                subp = csp[min_i]
                                del csp[min_i]
                                j = min_j
                                if min_t in [0,1]:
                                    if min_t == 0: j=j-1
                                    subp[-1][2], subp[0][0] = subp[-1][1], subp[0][1]
                                    subp = [ [subp[j][1], subp[j][1], subp[j][2]] ] + subp[j+1:] + subp[:j] + [ [subp[j][0], subp[j][1], subp[j][1]] ]
                                else:
                                    sp1,sp2,sp3 = csp_split(subp[j-1],subp[j],min_t)
                                    subp[-1][2], subp[0][0] = subp[-1][1], subp[0][1]
                                    subp = [ [ sp2[1], sp2[1],sp2[2] ] ] + [sp3] + subp[j+1:] + subp[:j-1] + [sp1] + [[ sp2[0], sp2[1],sp2[1] ]]
                                csp = [subp] + csp
                                if csp_subpath_ccw(csp[0]) :
                                    for i in range(len(csp)):
                                        n = []
                                        for j in csp[i]:
                                            n = [  [j[2][:],j[1][:],j[0][:]]  ] + n
                                        csp[i] = n[:]

                            d = cubicsuperpath.formatPath(csp)

                            d = re.sub(r'(?i)(m[^mz]+)',r'\1 Z ',d)
                            d = re.sub(r'(?i)\s*z\s*z\s*',r' Z ',d)
                            d = re.sub(r'(?i)\s*([A-Za-z])\s*',r' \1 ',d)

                        radius = self.unittouu(self.options.dotsize)*in_out
            # otp=path_modif()
                        nouvelle_courbe =inkex.etree.SubElement(path.getparent(), inkex.addNS('path','svg'),
                                {   inkex.addNS('label','inkscape'): 'gabarit: ' + self.options.dotsize + " " + self.options.stylegab,
                                    inkex.addNS('type','sodipodi'):	'inkscape:offset',
                                    inkex.addNS('radius','inkscape'):	str(radius),
                                    inkex.addNS('original','inkscape'):	d,
                                    'style': simplestyle.formatStyle({ 'stroke': '#00ff00', 'fill': 'none' ,'stroke-width': str(self.unittouu('1px'))}),
                                })
                    else:
                        inkex.errormsg("This extension need a path, not groups.")
            # otp.object_to_path(nouvelle_courbe)

if __name__ == '__main__':
    e = Dots()
    e.affect()


# vim: expandtab shiftwidth=4 tabstop=8 softtabstop=4 fileencoding=utf-8 textwidth=99
