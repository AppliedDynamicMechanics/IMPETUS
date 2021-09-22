// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha

#ifndef HEAD_DEF
#define HEAD_DEF
#include <math.h> 
#include <mpi.h> 

#define PI 3.14159265
#define DIM 3
typedef double real;
typedef double preal;
#define sqr(x) ((x)*(x))

typedef struct {
	unsigned int id;
	real m;
	preal x[DIM];
	preal x_original[DIM];
	preal v[DIM];
	preal F[DIM];
	preal nx[DIM];
	preal nv[DIM];
	preal na[DIM];
	int boundary_cross_count[DIM];
	bool islocal;
} Particle;

struct vec3{
	long double r[3];
};
struct TransformationMatrix {
	
	double H[3][3];
	double Hinv[3][3];
	TransformationMatrix(){
		for (int i = 0 ; i<3; i++){
			for (int j = 0 ; j < 3 ; j++ ) {
				H[i][j] = 0;
				Hinv[i][j] = 0;
			}
		}
		
	}
	void compute_inverse() {
		double DetH = 	H[0][0] * H[1][1] * H[2][2]+ H[0][1] * H[1][2] * H[2][0]+
						H[1][0] * H[2][1] * H[0][2]- H[0][2] * H[1][1] * H[2][0]-
						H[0][1] * H[1][0] * H[2][2]- H[1][2] * H[2][1] * H[0][0];

		Hinv[0][0] = (H[1][1] * H[2][2]- H[2][1] * H[1][2])/ DetH;
		Hinv[0][1] = (H[0][2] * H[2][1]- H[0][1] * H[2][2])/ DetH;
		Hinv[0][2] = (H[0][1] * H[1][2]- H[1][1] * H[0][2])/ DetH;
		Hinv[1][0] = (H[1][2] * H[2][0]- H[2][2] * H[1][0])/ DetH;
		Hinv[1][1] = (H[0][0] * H[2][2]- H[2][0] * H[0][2])/ DetH;
		Hinv[1][2] = (H[0][2] * H[1][0]- H[0][0] * H[1][2])/ DetH; 
		Hinv[2][0] = (H[1][0] * H[2][1]- H[1][1] * H[2][0])/ DetH;
		Hinv[2][1] = (H[0][1] * H[2][0]- H[0][0] * H[2][1])/ DetH;  
		Hinv[2][2] = (H[1][1] * H[0][0]- H[0][1] * H[1][0])/ DetH;
	}
	
	void VScaleCoordinate(preal x[DIM]) {
		preal x_T[DIM];
		x_T[0] = x[0] * Hinv[0][0] + x[1] * Hinv[1][0] + x[2] * Hinv[2][0];
		x_T[1] = x[0] * Hinv[0][1] + x[1] * Hinv[1][1] + x[2] * Hinv[2][1];
		x_T[2] = x[0] * Hinv[0][2] + x[1] * Hinv[1][2] + x[2] * Hinv[2][2];
		
		x[0] = x_T[0];
		x[1] = x_T[1];
		x[2] = x_T[2];
	}
	void VUnScaleCoordinate(preal x[DIM]) {
		preal x_T[DIM];
		x_T[0] = x[0] * H[0][0] + x[1] * H[1][0] + x[2] * H[2][0];
		x_T[1] = x[0] * H[0][1] + x[1] * H[1][1] + x[2] * H[2][1];
		x_T[2] = x[0] * H[0][2] + x[1] * H[1][2] + x[2] * H[2][2];
		x[0] = x_T[0];
		x[1] = x_T[1];
		x[2] = x_T[2];
		
	}
	void transform(preal x[DIM]){
		VUnScaleCoordinate (x);
	}
	void inverse_transform(preal x[DIM]){
		VScaleCoordinate (x);
	}
	void compute_inverse_and_inverse_transform(preal x[DIM]){
		compute_inverse();
		inverse_transform(x);
		//~ VScaleCoordinate (x);
	}
	
};
typedef struct ParticleList {
	Particle p;
	struct ParticleList *next;
} ParticleList;

typedef ParticleList* Cell;

#if 1==DIM
#define cell_index(ic,nc) ((ic)[0])
#elif 2==DIM
#define cell_index(ic,nc) ((ic)[0] + (nc)[0]*(ic)[1])
#elif 3==DIM
#define cell_index(ic,nc) ((ic)[0] + (nc)[0]*((ic)[1] + (nc)[1]*(ic)[2]))
#endif

#if 1==DIM
#define inverseindex(i,nc,ic) \
((ic)[0]=(i))
#elif 2==DIM
#define inverseindex(i,nc,ic) \
((ic)[0]=(i)%(nc)[0], (ic)[1]=(i)/(nc)[0])
#elif 3==DIM
#define inverseindex(i,nc,ic) \
((ic)[0]=(i)%(nc)[0], \
(ic)[1]=((i)/(nc)[0])%(nc)[1], \
(ic)[2]=((i)/(nc)[0])/(nc)[1])
#endif

#if 1==DIM
#define iterate(ic,minnc,maxnc) \
for ((ic)[0]=(minnc)[0]; (ic)[0]<(maxnc)[0]; (ic)[0]++)
#elif 2==DIM
#define iterate(ic,minnc,maxnc) \
for ((ic)[0]=(minnc)[0]; (ic)[0]<(maxnc)[0]; (ic)[0]++) \
for ((ic)[1]=(minnc)[1]; (ic)[1]<(maxnc)[1]; (ic)[1]++)
#elif 3==DIM
#define iterate(ic,minnc,maxnc) \
for ((ic)[0]=(minnc)[0]; (ic)[0]<(maxnc)[0]; (ic)[0]++) \
for ((ic)[1]=(minnc)[1]; (ic)[1]<(maxnc)[1]; (ic)[1]++) \
for ((ic)[2]=(minnc)[2]; (ic)[2]<(maxnc)[2]; (ic)[2]++)
#endif

inline void insertList(ParticleList **root_list, ParticleList *i) {
	i->next = *root_list;
	*root_list = i;
}

inline void deleteList(ParticleList **q) {
	*q = (*q)->next;
	 // (*q)->next points to element to be removed
}

inline int lengthList(ParticleList *i) {
	int counter = 0;
	while (NULL != i) {
		counter ++;
		i = i->next;
	}
	return counter;
}

inline preal v_dot(preal v1[DIM],preal v2[DIM]) {
	
	preal sum = 0;
	
	for (int d=0; d<DIM; d++) {
		sum += v1[d] * v2[d];
	}
	
	return sum;
}
inline double sum_mpi_double_all(MPI_Comm mpicomm, double local){
    double sum_global;

    MPI_Allreduce(&local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 
              mpicomm);
	return sum_global;
}
inline double sum_mpi_double(MPI_Comm mpicomm, double local, int receivingRank){
    double sum_global;
    MPI_Reduce(&local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, receivingRank, 
              mpicomm);
	return sum_global;
}
inline int sum_mpi_int(MPI_Comm mpicomm, int local, int receivingRank){
    int sum_global;
    MPI_Reduce(&local, &sum_global, 1, MPI_INT, MPI_SUM, receivingRank, 
              mpicomm);
	return sum_global;
}

#endif
