from __future__ import division
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
from scikits.sparse.cholmod import cholesky
from scipy.sparse import csc_matrix, diags, linalg as sla
import numpy as np
import vtk
import timeit

# filename = "SourceVtp/largeMesh.vtp"
# outputfilename = "largeRes.vtp"
filename = "SourceVtp/simple.vtp"
outputfilename = "largeRes.vtp"
samplenum = 100

def check_correlation(X, npNodes):
	ptsidx = np.arange(11) * 12 # choose points on the diagonal
	ptsidx[1] = 2 # adjust with the real point lable

	# ptsidx = np.arange(101) * 2 # choose points along the first row
	# ptsidx[1] = 1 # adjust with the real point lable
	corX = np.corrcoef(X)
	distance = np.sqrt(np.sum((npNodes[ptsidx] - npNodes[0]) ** 2, axis=1))
	plt.plot(distance, corX[0, ptsidx])
	plt.ylabel('Correlation')
	plt.xlabel('Distance')
	plt.show()

def loc(indptr, indices, i, j):
	return indptr[i] + np.where(indices[indptr[i]:indptr[i+1]]==j)[0]

def main():

	kappa = 0.5656

	start_time = timeit.default_timer()

	# Read triangulation from file.
	print 'Reading File...'
	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(filename)
	reader.Update()
	polyDataModel = reader.GetOutput()

	totalNodes = polyDataModel.GetNumberOfPoints()
	vtkNodes = polyDataModel.GetPoints().GetData()
	npNodes = vtk_to_numpy(vtkNodes)
	# print(vtkPointData)
	# print(npNodes)

	print 'Total nodes: ', totalNodes
	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Building Topology...'
	totalElms = polyDataModel.GetNumberOfCells()
	# Get cells from source file.
	npElms = np.zeros((totalElms, 3), dtype=int)
	npEdges = np.zeros((totalElms, 3, 3))
	npAreas = np.zeros(totalElms)
	for icell in xrange(totalElms):
		cell = polyDataModel.GetCell(icell)
		numpts = cell.GetNumberOfPoints()
		for ipt in xrange(numpts):
			npElms[icell, ipt] = cell.GetPointId(ipt)
		npEdges[icell, 0] = npNodes[npElms[icell, 2]] - npNodes[npElms[icell, 1]]
		npEdges[icell, 1] = npNodes[npElms[icell, 0]] - npNodes[npElms[icell, 2]]
		npEdges[icell, 2] = npNodes[npElms[icell, 1]] - npNodes[npElms[icell, 0]]
		npAreas[icell] = cell.ComputeArea()
		# for iedge in xrange(numedges):
		# 	edge = cell.GetEdge(iedge)
		# 	print edge
		# 	return
	# print(npElms)
	# print(npEdges)
	# print(npAreas)

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Creating the sparse matrix...'
	# Create sparse data structure.
	sparseInfo = [[] for _ in xrange(totalNodes)]
	for icell in xrange(totalElms):
		for inode in npElms[icell]:
			# [sparseInfo[inode].extend([pt]) for pt in npElms[icell] if pt not in sparseInfo[inode]]
			sparseInfo[inode].extend(npElms[icell])
	sparseInfo = np.array(sparseInfo)
	for knodes in xrange(totalNodes):
		sparseInfo[knodes] = np.unique(sparseInfo[knodes])
	# print(sparseInfo)

	# Generate the sparse matrix.
	indptr = [0]
	indices = []
	for inode in xrange(totalNodes):
		indices.extend(sparseInfo[inode])
		indptr.extend([len(indices)])
	rawC = np.zeros(len(indices))
	rawG = np.zeros(len(indices))
	rawCTuta = np.zeros(totalNodes)

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Assembling global matrix...'
	# Generate C and G matrix.
	cm = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]])
	dcm = np.array([1.0, 1.0, 1.0])
	for icell in xrange(totalElms):
		# Compute local matrix first.
		localc = cm * npAreas[icell] / 12.0
		localdc = dcm * npAreas[icell] / 3.0
		localg = np.dot(npEdges[icell], npEdges[icell].transpose()) / (4 * npAreas[icell])
		# Assembly to the glabal matrix.
		for i in xrange(3):
			rawCTuta[npElms[icell, i]] += localdc[i]
			for j in xrange(3):
				rawindex = loc(indptr, indices, npElms[icell, i], npElms[icell, j])
				rawC[rawindex] += localc[i, j]
				rawG[rawindex] += localg[i, j]

	C = csc_matrix((rawC, np.array(indices), np.array(indptr)), shape=(totalNodes, totalNodes))
	G = csc_matrix((rawG, np.array(indices), np.array(indptr)), shape=(totalNodes, totalNodes))

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Creating inverse C...'
	invCTuta = diags([1.0 / rawCTuta], [0], shape=(totalNodes, totalNodes))

	# print 'Computating C Inverse...'
	# factorC = cholesky(C)
	# invC = factorC.inv()

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	# Compute Q matrix according to C and G.
	print 'Computing K...'
	K = (kappa**2)*C + G

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Computing of Q...'
	Q1 = K
	# Q = (K.dot(invC)).dot(K)
	# Q = (K.dot(invCTuta)).dot(K) # Q2
	Q = (((K.dot(invCTuta)).dot(Q1)).dot(invCTuta)).dot(K)
	# Q = (((K.dot(invC)).dot(Q1)).dot(invC)).dot(K)

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Cholesky factor of Q...'
	# Decomposition.
	factorQ = cholesky(Q)
	# L = factorQ.L()
	# print(factorQ.L())
	# lu = sla.splu(Q)
	# print(lu.L)

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Generating samples...'
	# Generate normal distrib random nums & combine.
	Z = np.random.normal(size=(totalNodes, samplenum))

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'Solving upper triangular syms...'
	X = factorQ.solve_Lt(Z)

	print timeit.default_timer() - start_time
	start_time = timeit.default_timer()

	print 'std: ', np.std(X[1:4], axis=1)

	# Store back the random field.
	print 'Exporting data...'
	vtkPointData = polyDataModel.GetPointData()
	for itrade in xrange(100): # X.shape[1] # !not exporting all data to save time
		scaler = numpy_to_vtk(np.ascontiguousarray(X[:,itrade]))
		scaler.SetName('RandomField ' + str(itrade+1))
		vtkPointData.AddArray(scaler)

	writer = vtk.vtkXMLPolyDataWriter()
	writer.SetInputData(polyDataModel)
	writer.SetFileName(outputfilename)
	writer.Write()

	print timeit.default_timer() - start_time

	# print 'Ploting...'
	# check_correlation(X, npNodes)

	# ----------------------------------------
	# np.savetxt('large-nodes.out', npNodes)
	# np.savetxt('large-X-invCTuta.out', X[:,:100])
	# # Save the sparse info into file
	# # np.savetxt('simple-sparseInfo.out', sparseInfo, fmt='%s')
	# with open("simple-sparseInfo.out","w") as f:
	# 	f.write("\n".join(" ".join(map(str, x)) for x in sparseInfo))
	# ----------------------------------------

if __name__ == '__main__':
	main()