# Jain, Kritika
# 1001-093-381
# 2017-09-24
# Assignment_02_03

import numpy
import math

class cl_world:
    def __init__(self, objects=[], canvases=[]):
        self.objects = objects
        self.canvases = canvases
        # variables are initialized
        self.vertexes = []
        self.faces = []
        self.window = []
        self.viewport = []

    def add_canvas(self, canvas):
        self.canvases.append(canvas)
        canvas.world = self


    def map_to_viewport(self,pixels,viewport,window):
        for point in pixels:
            for i in range(0,2):
                point[i] = (point[i]*(viewport[i+2] - viewport[i])/(window[i+2]-window[i]) + (viewport[i+2]+viewport[i])/2) #/(viewport[i+2] - viewport[i])#)*/(viewport[i+2] - viewport[i])#*  #(viewport[i+2] - viewport[i])# + viewport[i]  # this is how maping to wiewport is performed
        return pixels


    def map_XY(self, vertexes):
        result = []
        for point in vertexes:
            result.append([point[0],point[1]])  # proecting to xy
        result = result
        return result


    def map_to_window(self, pixels, window, canvas):
        for point in pixels:
            point[0] = int(int(canvas.cget("width"))*((point[0]*(window[2]-window[0]) +(window[2]+window[0])/2)/(window[2]-window[0])))
            # next line accounts for the direction of y
            point[1] = int(canvas.cget("height"))-int(int(canvas.cget("height")) * (((point[1])*(window[3]-window[1])+(window[3]+window[1])/2)/(window[3]-window[1]))) ##cget("height"))-
        return pixels

    def prepare_port_view(self, polynom , canvas, window):
        polynom1 = [0,0,0,0]
        polynom1[0] = int(int(canvas.cget("width"))* polynom[0] )
        polynom1[1] =int(canvas.cget("height"))-  int(int(canvas.cget("height"))*polynom[1] )
        polynom1[2] = int(int(canvas.cget("width")) * polynom[2] )
        polynom1[3] =int(canvas.cget("height"))- int(int(canvas.cget("height")) * polynom[3] )
        result = [polynom1[0],polynom1[1],polynom1[2],polynom1[1],polynom1[2],polynom1[3],polynom1[0],polynom1[3]]
        return result


    def create_graphic_objects(self, canvas):
        pixels =self.map_XY(self.vertexes)
        pixels = self.map_to_viewport(pixels, self.viewport, self.window)
        pixels=self.map_to_window(pixels, self.window, canvas)

        viewport_frame = self.prepare_port_view(self.viewport, canvas,self.window)
        self.objects = []
        self.objects.append(canvas.create_polygon(viewport_frame, fill="yellow", outline='black'))
        for face in self.faces:
            # this is how we add polys to canvas
            self.objects.append(canvas.create_polygon(pixels[face[0]-1][0],pixels[face[0]-1][1],
                                                      pixels[face[1]-1][0],pixels[face[1]-1][1],
                                                      pixels[face[2]-1][0],pixels[face[2]-1][1],
                                                      fill="red",outline='black'))


    def redisplay(self, canvas, ):
        pixels = self.map_XY(self.vertexes)
        pixels = self.map_to_viewport(pixels, self.viewport,self.window)
        pixels = self.map_to_window(pixels, self.window, canvas)
        for face in self.faces:
            # in here only the coordinates of polys are changed, so they are redrawn
            new_coords =[pixels[face[0]-1][0],pixels[face[0]-1][1],
                                                      pixels[face[1]-1][0],pixels[face[1]-1][1],
                                                      pixels[face[2]-1][0],pixels[face[2]-1][1]]
            canvas.coords(self.objects[self.faces.index(face)+1], new_coords)

    def rotate(self, angle, start_point, end_point, vertexes):
        def rotation_matrix(vector, theta):  # rotation matrix constructed here
            vector = numpy.asarray(vector)  
            vector = vector / math.sqrt(numpy.dot(vector, vector))
            a = math.cos(theta / 2.0)
            b, c, d = -vector * math.sin(theta / 2.0)
            aa, bb, cc, dd = a * a, b * b, c * c, d * d
            bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
            return numpy.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                                 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                                 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

        vector = [0, 0, 0]
        for i in range(3):
            vector[i] = float(end_point[i]) - float(start_point[i]) # a vector for rotation constructed
        vectro_length = (vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2) ** (0.5) # vector length
        for i in range(3):  # vector normalized
            vector[i] = vector[i] / (vectro_length)
        vertexes1 =[]
        angle = angle * math.pi/180
        for vertex in vertexes:  # rotation is made as multiplication of matrices
            vertex1 = numpy.dot(rotation_matrix(vector, angle),vertex)
            vertexes1.append(vertex1)
        return vertexes1

    def translate_vector(self, translate, vertexes):
        for vertex in vertexes:
            for i in range(3):  # translation is a simple shift of coordinates
                vertex[i] = vertex[i] + translate[i]
        return vertexes

    def scale(self, vector, scale, vertexes): # scale to point is a combination of a scale over origin
        for vertex in vertexes:  # and shift
            for i in range(3):  # first we subtract the point
                vertex[i] = vertex[i] - vector[i]
        scale_matrix = [[scale[0], 0, 0], [0, scale[1], 0], [0, 0, scale[2]]]  # scaling matrix
        vertexes1 = []
        for vertex in vertexes:
            vertex1 = numpy.dot(scale_matrix, vertex)
            vertexes1.append(vertex1)
        for vertex in vertexes1:
            for i in range(3):  # and then we add it (after scale)
                vertex[i] = vertex[i] + vector[i]
        return vertexes1

