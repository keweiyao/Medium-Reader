# cython: c_string_type=str, c_string_encoding=ascii
from libc.math cimport *
from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libcpp cimport bool
import h5py
import numpy as np
cimport numpy as np

cdef inline double finterp(double *** c, double rx, double ry, double rz):
	cdef double result = 0.0
	cdef size_t i, j, k
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
	cdef object _f, _step_keys, info_keys, static_property, _tabs
	cdef str _mode
	cdef public size_t _Nx, _Ny, _step_key_index
	cdef public double _dx, _dy, _xmin, _ymin, _xmax, _ymax, _tstart, _dt, _tnow
	cdef double T_static
	cdef double *** interpcube 
	cdef bool status
	
	def __cinit__(self, medium_flags):
		self._mode = medium_flags['type']
		
		if self._mode == 'static':
			print "works in static meidum mode!"
			print "Medium property can be specified step by step"
			self._Nx = 0
			self._Ny = 0
			self._dx = 0.
			self._dy = 0.
			self._tstart = 0.0
			self._dt = medium_flags['static_dt']
			self._xmin = 0.
			self._xmax = 0.
			self._ymin = 0.
			self._ymax = 0.
			self._tnow = self._tstart - self._dt
			self.status = True
		elif self._mode == 'dynamic':
			hydrofilename = medium_flags['hydrofile']
			if hydrofilename == None:
				raise ValueError("Need hydro history file")
			else:
				self._f = h5py.File(hydrofilename, 'r')

			self._step_keys = list(self._f['Event'].keys())
			self._step_key_index = 0
			self.status = True
			self.info_keys = self._f['Event'][self._step_keys[0]].keys() 
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
			self._tnow = self._tstart - self._dt
			
			self.interpcube = <double ***> malloc(2*sizeof(double**))
			for i in range(2):
				self.interpcube[i] = <double **> malloc(2*sizeof(double*))
				for j in range(2):
					self.interpcube[i][j] = <double *> malloc(2*sizeof(double))
		else:
			raise ValueError("Medium mode not implemented.")
	
	cpdef init_tau(self):
		return self._tstart
	cpdef hydro_status(self):
		return self.status
	cpdef dtau(self):
		return self._dt
	cpdef boundary(self):
		return self._xmin, self._xmax, self._ymin, self._ymax
	
	cdef frame_inc_unpack(self):
		self._tabs = {}
		cdef vector[vector[vector[double]]] buff
		buff.resize(2)
		for key2 in self.info_keys:
			for i in range(2):
				key1 = self._step_keys[self._step_key_index+i]
				buff[i] = self._f['Event'][key1][key2].value
			self._tabs.update({key2:buff})
		self._step_key_index += 1
		if self._step_key_index == len(self._step_keys) - 1:
			self.status = False

	cpdef load_next(self, StaticPropertyDictionary=None):
		if self._mode == "dynamic":
			self.frame_inc_unpack()
		elif self._mode == 'static':
			if StaticPropertyDictionary == None:
				raise ValueError("Need to provide a static property at this step")
			else:			
				self.static_property = StaticPropertyDictionary
		else:
			raise ValueError("Medium mode not implemented.")
		self._tnow += self._dt
	
	cpdef get_current_frame(self, key):
		return np.array(self._tabs[key][0])
	
	cpdef interpF(self, double tau, xvec, keys):
		cdef double rt, nx, ny, rx, ry, gamma, buff, vz
		cdef int ix, iy, i, j, k
		if self._mode == "static":
			return [self.static_property[key] for key in keys]
		if self._mode == "dynamic":
			if xvec[1] < self._xmin or xvec[1] > self._xmax \
				or xvec[2] < self._ymin or xvec[2] > self._ymax:
				return [0. for key in keys]
			else:
				result = []
				rt = (tau - self._tnow)/self._dt
				nx = (xvec[1] - self._xmin)/self._dx
				ny = (xvec[2] - self._ymin)/self._dy
				ix = <int>floor(nx)
				rx = nx - ix
				iy = <int>floor(ny)
				ry = ny - iy
				vz = xvec[3]/xvec[0]
				for key in keys:
					if key == 'Vz':
						result.append(vz)
					else:
						for k in range(2):
							for i in range(2):
								for j in range(2):
									self.interpcube[k][i][j] = self._tabs[key][k][ix+i][iy+j]
						buff = finterp(self.interpcube, rt, rx, ry)
						if key == 'Vx' or key == 'Vy':
							gamma = 1.0/sqrt(1.0-vz*vz)
							buff /= gamma
						result.append(buff)
				return result
		


		
		

	
