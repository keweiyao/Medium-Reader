from medium cimport *
from libc.math cimport *
from libc.stdio cimport *
from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np
import h5py
from cpython cimport array
import array

cdef inline finterp(double *** c, double rx, double ry, double rz):
	cdef double result = 0.0
	cdef double vx[2]
	vx[0] = 1. - rx; vx[1] = rx
	cdef double vy[2]
	vy[0] = 1. - ry; vy[1] = ry
	cdef double vz[2]
	vz[0] = 1. - rz; vz[1] = rz
	for i in range(2):
		for j in range(2):
			for k in range(2):
				result += c[i][j][k]*vx[i]*vy[j]*vz[k]
	return result

cdef class Medium:
	cdef object _f, _keys, _tabs0, _tabs1
	cdef public size_t _Nx, _Ny, _key_index
	cdef public double _dx, _dy, _xmin, _ymin, _xmax, _ymax, _tstart, _dt, _tnow
	cdef double *** interpcube
	
	def __cinit__(self, filename):
		self._f = filename
		try:
			self._f = h5py.File(filename, 'r')
		except IOError:
			print "Open hydrofile %s failed."%filename
		self._keys = self._f['Event'].keys()
		self._Nx = self._f['Event'].attrs['XH'] - self._f['Event'].attrs['XL'] + 1
		self._Ny = self._f['Event'].attrs['YH'] - self._f['Event'].attrs['YL'] + 1
		self._dx = self._f['Event'].attrs['DX']
		self._dy = self._f['Event'].attrs['DY']
		self._tstart = self._f['Event'].attrs['Tau0']
		self._dt = self._f['Event'].attrs['dTau']
		self._xmin = self._f['Event'].attrs['XL']*self._dx
		self._xmax = self._f['Event'].attrs['XH']*self._dx
		self._ymin = self._f['Event'].attrs['YL']*self._dy
		self._ymax = self._f['Event'].attrs['YH']*self._dy
		self._key_index = 0
		self._tnow = self._tstart - self._dt

		self.interpcube = <double ***> malloc(2*sizeof(double**))
		for i in range(2):
			self.interpcube[i] = <double **> malloc(2*sizeof(double*))
			for j in range(2):
				self.interpcube[i][j] = <double *> malloc(2*sizeof(double))

	cpdef unpack_frame(self, index):
		key = self._keys[index]
		T  = self._f['Event'][key]['Temp'].value
		Vx = self._f['Event'][key]['Vx'].value
		Vy = self._f['Event'][key]['Vy'].value
		e = self._f['Event'][key]['e'].value
		p = self._f['Event'][key]['P'].value
		pi00 = self._f['Event'][key]['Pi00'].value
		pi01 = self._f['Event'][key]['Pi01'].value
		pi02 = self._f['Event'][key]['Pi02'].value
		pi03 = self._f['Event'][key]['Pi03'].value
		pi11 = self._f['Event'][key]['Pi11'].value
		pi12 = self._f['Event'][key]['Pi12'].value
		pi13 = self._f['Event'][key]['Pi13'].value
		pi22 = self._f['Event'][key]['Pi22'].value
		pi23 = self._f['Event'][key]['Pi23'].value
		pi33 = self._f['Event'][key]['Pi33'].value
		return T, Vx, Vy, e, p, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33

	cpdef load_next(self):
		if self._key_index == 0:
			T, Vx, Vy, e, p, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33 \
											 = self.unpack_frame(self._key_index)
			self._tabs0 = {'Temp': T, 'Vx': Vx, 'Vy': Vy, 'e': e, 'p': p, 
						   'pi00': pi00, 'pi01': pi01, 'pi02': pi02, 'pi03': pi03, 
						   'pi11': pi11, 'pi12': pi12, 'pi13': pi13, 'pi22': pi22, 
						   'pi23': pi23, 'pi33': pi33}
		else:
			self._tabs0 = self._tabs1
		self._key_index += 1
		self._tnow += self._dt
		T, Vx, Vy, e, p, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33 \
											 = self.unpack_frame(self._key_index)
		self._tabs1 = {'Temp': T, 'Vx': Vx, 'Vy': Vy, 'e': e, 'p': p, 
					   'pi00': pi00, 'pi01': pi01, 'pi02': pi02, 'pi03': pi03, 
					   'pi11': pi11, 'pi12': pi12, 'pi13': pi13, 'pi22': pi22, 						   'pi23': pi23, 'pi33': pi33}
	
	
	cpdef double interpF(self, t, xvec, key):
		if xvec[0] < self._xmin or xvec[0] > self._xmax \
			or xvec[1] < self._ymin or xvec[1] > self._ymax:
			return 0
		cdef rt = (t - self._tnow)/self._dt
		cdef double nx = (xvec[0] - self._xmin)/self._dx
		cdef double ny = (xvec[1] - self._ymin)/self._dy
		cdef int ix = <int>floor(nx)
		cdef rx = nx - ix
		cdef int iy = <int>floor(ny)
		cdef ry = ny - iy
		cdef int i, j
		for i in range(2):
			for j in range(2):
				self.interpcube[0][i][j] = self._tabs0[key][ix+i, iy+i]
				self.interpcube[1][i][j] = self._tabs1[key][ix+i, iy+i]

		return finterp(self.interpcube, rt, rx, ry)
		


		
		

	
