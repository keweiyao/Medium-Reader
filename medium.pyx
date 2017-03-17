# cython: c_string_type=str, c_string_encoding=ascii
from medium cimport *
from libc.math cimport *
from libc.stdio cimport *
from libc.stdlib cimport malloc, free
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
	cdef object _f, _step_keys, _tabs0, _tabs1, _mode
	cdef public size_t _Nx, _Ny, _step_key_index
	cdef public double _dx, _dy, _xmin, _ymin, _xmax, _ymax, _tstart, _dt, _tnow
	cdef double T_static
	cdef double *** interpcube
	
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
		elif self._mode == 'dynamic':
			hydrofilename = medium_flags['hydrofile']
			if hydrofilename == None:
				raise ValueError("Please provide the hydro histroy file in dynamic meidum model.")
			try:
				self._f = h5py.File(hydrofilename, 'r')
			except IOError:
				print "Open hydrofile %s failed."%hydrofilename
			self._step_keys = self._f['Event'].keys()
			self._step_key_index = 0
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
			raise ValueError("Medium mode %s not implemented."%self._mode)
	
	cpdef init_tau(self):
		return self._tstart
	cpdef dtau(self):
		return self._dt
	cpdef boundary(self):
		return self._xmin, self._xmax, self._ymin, self._ymax
	
	cpdef unpack_frame(self, step_index, staticT = None):
		key = self._step_keys[step_index]
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
			
	cpdef load_next(self, StaticPropertyDictionary=None):
		status = True
		if self._mode == "dynamic":
			if self._step_key_index == 0:
				T, Vx, Vy, e, p, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33 \
												 = self.unpack_frame(self._step_key_index)
				self._tabs0 = {'Temp': T, 'Vx': Vx, 'Vy': Vy, 'e': e, 'p': p, 
							   'pi00': pi00, 'pi01': pi01, 'pi02': pi02, 'pi03': pi03, 
							   'pi11': pi11, 'pi12': pi12, 'pi13': pi13, 'pi22': pi22, 
							   'pi23': pi23, 'pi33': pi33}
			else:
				self._tabs0 = self._tabs1
			self._step_key_index += 1
			T, Vx, Vy, e, p, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33 \
												 = self.unpack_frame(self._step_key_index)
			self._tabs1 = {'Temp': T, 'Vx': Vx, 'Vy': Vy, 'e': e, 'p': p, 
						   'pi00': pi00, 'pi01': pi01, 'pi02': pi02, 'pi03': pi03, 
						   'pi11': pi11, 'pi12': pi12, 'pi13': pi13, 'pi22': pi22,								'pi23': pi23, 'pi33': pi33}
			if self._step_key_index == len(self._step_keys) - 1:
				status = False
		if self._mode == 'static':
			if StaticPropertyDictionary == None:
				raise ValueError("Need to provide a static property at this step")
			self._tabs0 = StaticPropertyDictionary
		self._tnow += self._dt
		return status
	
	cpdef get_current_frame(self, key):
		return self._tabs0[key]
	
	cpdef interpF(self, tau, xvec, keys):
		cdef double rt, nx, ny, rx, ry, gamma, buff
		cdef int ix, iy, i, j
		if self._mode == "static":
			return [self._tabs0[key] for key in keys]
		if self._mode == "dynamic":
			if xvec[0] < self._xmin or xvec[0] > self._xmax \
				or xvec[1] < self._ymin or xvec[1] > self._ymax:
				return [0.0]*len(keys)
			rt = (tau - self._tnow)/self._dt
			nx = (xvec[0] - self._xmin)/self._dx
			ny = (xvec[1] - self._ymin)/self._dy
			ix = <int>floor(nx)
			rx = nx - ix
			iy = <int>floor(ny)
			ry = ny - iy
			result = []
			vz = xvec[2]/xvec[3]
			gamma = 1.0/sqrt(1.0-vz*vz)
			for key in keys:
				if key == 'Vz':
					result.append(vz)
				else:
					for i in range(2):
						for j in range(2):
							self.interpcube[0][i][j] = self._tabs0[key][ix+i, iy+i]
							self.interpcube[1][i][j] = self._tabs1[key][ix+i, iy+i]
					buff = finterp(self.interpcube, rt, rx, ry)
					#Note that this vx and vy are at mid-rapidity and need to be boosted to obtain the solution at forward and backward rapidity
					if key == 'Vx' or keys == 'Vy':
						buff /= gamma
					result.append(buff)
			return result
		


		
		

	
