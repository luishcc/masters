## =================================================================== ##
#  this is file InOut.py, created at 12-Jun-2013                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##

import numpy as np

class InOut:
 def __init__(_self,_X,_Y,_IEN,_numVerts,_numElems,_scalar , _scalar2=None, _vet1=None, _vet2=None):
  _self.X = _X
  _self.Y = _Y
  _self.IEN = _IEN
  _self.numVerts = _numVerts
  _self.numElems = _numElems
  _self.scalar = _scalar
  _self.scalar2 = _scalar2
  _self.vet1 = _vet1
  _self.vet2 = _vet2


 def saveVTK(_self,_dir,_file,_iter=None):

  if _iter is None:
   vtkFile = open(_dir + '/' + _file + '.vtk', 'w')
   _self.vtkHeader(vtkFile)
  else:
   vtkFile = open(_dir + '/' + _file + '-' + str(_iter) + '.vtk', 'w')
   _self.vtkHeader(vtkFile,_iter)

  _self.vtkCoords(vtkFile)
  _self.vtkCellArray(vtkFile)
  _self.vtkCellType(vtkFile)
  _self.vtkScalarScalarHeader(vtkFile)

  if _self.scalar is not None:
   _self.vtkScalar(vtkFile,"Omega",_self.scalar);

  if _self.scalar2 is not None:
   _self.vtkScalar(vtkFile,"Psi",_self.scalar2);

  if _self.vet1 is not None:
   _self.vtkVector(vtkFile,"velocidade",_self.vet1, _self.vet2);

  vtkFile.close()



 def vtkHeader(_self,_file,_iter=None):
  _file.write( "# vtk DataFile Version 1.0\n" )
  _file.write( "2D Simulation C++\n" )
  _file.write( "ASCII\n" )
  _file.write( "DATASET UNSTRUCTURED_GRID\n" )

  _file.write( "FIELD FieldData 1\n" )
  _file.write( "NODES 1 2 int\n" )
  _file.write( str(_self.numVerts) + " " + \
               str(_self.numElems) + "\n" )

  _file.write( "\n" )

 def vtkCoords(_self,_file):
  _file.write( "POINTS " + str(_self.numVerts) + " double\n" )
  for i in range(0,_self.numVerts):
   _file.write( str(_self.X[i]) + " " + \
                str(_self.Y[i]) + " 0.0\n" )

  _file.write( "\n" )

 def vtkCellArray(_self,_file):
  _file.write( "CELLS " + str(_self.numElems) \
                        + " " + str(4*_self.numElems) + "\n" )
  for i in range(0,_self.numElems):
   _file.write( "3 " + str(_self.IEN[i][0]) + " " + \
                       str(_self.IEN[i][1]) + " " + \
                       str(_self.IEN[i][2]) + "\n" )

  _file.write( "\n" )

 def vtkCellType(_self,_file):
  _file.write( "CELL_TYPES " + str(_self.numElems) + "\n" )
  for i in range(0,_self.numElems):
   _file.write( "5 ")

  _file.write( "\n" )
  _file.write( "\n" )
   
 def vtkScalarHeader(_self,_file):
  _file.write( "POINT_DATA " + str(_self.numVerts) + "\n" )
  _file.write( "\n" )

 def vtkScalarScalarHeader(_self,_file):
  _file.write( "POINT_DATA " + str(_self.numVerts) + "\n" )

 def vtkScalar(_self,_file,_name,_scalar):
  _file.write( "SCALARS " + _name + " double\n" )
  _file.write( "LOOKUP_TABLE default\n" )

  for i in range(0,_self.numVerts):
   _file.write( str(_scalar.item(i)) + "\n" )

  _file.write( "\n" )

 def vtkVector(_self,_file,_name,_vec1,_vec2):
  _file.write( "VECTORS " + _name + " double\n" )

  for i in range(0,_self.numVerts):
   _file.write( str(_vec1.item(i)) + " " + \
                str(_vec2.item(i)) + " 0.0\n" )

  _file.write( "\n" )



     
