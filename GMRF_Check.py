from __future__ import division
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import matplotlib.pyplot as plt
import numpy as np
import vtk
from shortestpath import Dijkstra

# filename = "SourceVtp/largeMesh.vtp"
# outputfilename = "largeMesh-100.vtp"

nodesfile = "simple-res-Q3/simple-nodes.out"
xfile = "simple-res-Q3/simple-X-invC.out" # invCTuta
Gfile = "simple-res-Q3/simple-sparseInfo.out"

# nodesfile = "large-res-Q3/large-nodes.out"
# xfile = "large-res-Q3/large-X-invCTuta.out"
# Gfile = "large-res-Q3/large-sparseInfo.out"

point_sample_num = 100

if __name__ == '__main__':

	# reader = vtk.vtkXMLPolyDataReader()
	# reader.SetFileName(outputfilename)
	# reader.Update()
	# polyDataModel = reader.GetOutput()

	# totalNodes = polyDataModel.GetNumberOfPoints() # total nodes number
	# vtkNodes = polyDataModel.GetPoints().GetData()
	# npNodes = vtk_to_numpy(vtkNodes) # nodes matrix

	# vtkPointDataScalars = polyDataModel.GetPointData().GetScalars('RandomField 1')
	# X = np.zeros(totalNodes)
	# for i in xrange(totalNodes):
	# 	X[i] = vtkPointDataScalars.GetValue(i)


	# # calculating the distance between the 1st node and all others
	# distances = np.linalg.norm(npNodes-npNodes[0,:], axis=1)
	# print len(distances)

	# # calculating the covariance matrix
	# corX = np.cov(X)

	# plt.plot(distances, corX[0,:])
	# plt.show()

	nodes = np.loadtxt(nodesfile)
	X = np.loadtxt(xfile)
	print X.shape[0]

	# # -------------- plot along a line --------------
	# # find points along the lower boundary
	# idx = np.where(nodes[:,1] == np.min(nodes[:,1]))

	# cnodes = nodes[idx]
	# cX = X[idx]
	
	# # distances
	# distances = np.linalg.norm(cnodes-cnodes[0], axis=1)
	# # covariances
	# # cov = np.cov(cX)
	# cor = np.corrcoef(cX)

	# # plot
	# plt.plot(distances, cor[0])
	# plt.show()

	# # -------------- plot randomly --------------
	# # generate sample indices
	# idx = np.random.choice(X.shape[0], point_sample_num)
	# cnodes = nodes[idx]
	# cX = X[idx]

	# # calculating the distance between the 1st node and all others
	# distances = np.linalg.norm(cnodes-cnodes[0], axis=1)
	# sortIndices = np.argsort(distances)

	# # calculating the covariance matrix
	# cor = np.corrcoef(cX)
	# # cov = np.cov(cX)

	# plt.plot(distances[sortIndices], cor[0,sortIndices])
	# plt.show()

	# -------------- plot randomly using distances along elements --------------
	G = [np.array(map(int, line.split())) for line in open(Gfile)]

	# # find a start and end point
	# idx = np.random.choice(X.shape[0], point_sample_num)
	# idx = np.sort(idx)

	# # calculate the path length
	# D,P = Dijkstra(nodes,G,idx[0],idx[-1])

	# # calculate distances from start to all the other points
	# distances = np.zeros(point_sample_num)
	# for i in xrange(len(idx)):
	# 	distances[i] = D[idx[i]]

	# sortIndices = np.argsort(distances)
	# # calculate correlatives
	# cor = np.corrcoef(X[idx])

	# plt.plot(distances[sortIndices], cor[0,sortIndices])
	# plt.show()

	# # all points with distances of grid
	# D,P = Dijkstra(nodes,G,0)
	# distances = np.zeros(len(nodes))
	# for i in xrange(len(nodes)):
	# 	distances[i] = D[i]
	# sortIndices = np.argsort(distances)
	# cor = np.corrcoef(X)
	# plt.plot(distances[sortIndices], cor[0,sortIndices])
	# plt.show()


