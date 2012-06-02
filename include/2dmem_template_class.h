#ifndef __2DMEM_TEMPLATE_CLASS_H__
#define __2DMEM_TEMPLATE_CLASS_H__

/**********************************************
* Add the operator for directly matrice operation *
*
*                                             *
*                 Written by Pei-Chi, Hsiao*
*                 Released 2009/06/08	14:48	*
**********************************************/

/**********************************************
* Add m_bFree for directly assigning a memory *
* created from other part.                    *
* The memory should be freed from other part. *
*                                             *
*                 Written by Kai-Yueh, Chang  *
*                 Released 2008/11/25 11:53   *
**********************************************/

/**********************************************
* Add the operator[] for directly accessing   *
*                                             *
*                 Written by Kai-Yueh, Chang  *
*                 Released 2006/04/19 07:40   *
**********************************************/

/**********************************************
* For the problem of the function reveal in   *
* VC++ 6.0, the functions are moved into      *
* class bady.                                 *
*                                             *
* For acceleration, some functions are        *
* changed to be inline function.              *
*                                             *
*                 Written by Kai-Yueh, Chang  *
*                 Released 2004/04/06 20:02   *
**********************************************/

/**********************************************
* This template class can allocate 2D memory. *
*                                             *
* Copyleft. You can copy it and modify it     *
* suitably for your purposes. And, if you     *
* have any question, don't ask or e-mail to   *
* me, solve it by yourself or consult with    *
* others.                                     *
*                                             *
*                 Written by Kai-Yueh, Chang  *
*                 Released 2003/07/14 00:04   *
**********************************************/

#include <stdlib.h>
#include <math.h>

template<class T>
class TC2DMem{
public:
	int m_row;
	int m_col;
	T** m_ppMem;
private:
//public:
	bool m_bFree;

public:
	TC2DMem(){
		m_row = 0;
		m_col = 0;
		m_bFree = 0;
		m_ppMem = NULL;
	}
	~TC2DMem(){
		Free2DMem();
	}

	T** Alloc2DMem(int nRow, int nCol){
		T** ppMem;
		T* pMem;

		T** ppCur;
		T* pCur;

		int tmp;

		Free2DMem();
		if(nRow<1)
			return NULL;
		if(nCol<1)
			return NULL;
		m_row = nRow;
		m_col = nCol;

		ppMem = (T**)malloc(m_row*sizeof(T*));
		pMem = new T[m_row*m_col];
		ppCur = ppMem;
		pCur = pMem;

		for(tmp=0; tmp<m_row; tmp++){
			*ppCur = pCur;
			pCur += m_col;
			ppCur++;
		}

		m_ppMem = ppMem;

		m_bFree = 1;
		return m_ppMem;
	}
	T** AssignMemory(int nRow, int nCol, T* pMem){
		T** ppMem;

		T** ppCur;
		T* pCur;

		int tmp;

		Free2DMem();
		if(nRow<1)
			return NULL;
		if(nCol<1)
			return NULL;
		m_row = nRow;
		m_col = nCol;

		ppMem = (T**)malloc(m_row*sizeof(T*));
		ppCur = ppMem;
		pCur = pMem;

		for(tmp=0; tmp<m_row; tmp++){
			*ppCur = pCur;
			pCur += m_col;
			ppCur++;
		}

		m_ppMem = ppMem;

		return m_ppMem;
	}
	void Free2DMem(){
		if(NULL!=(m_ppMem)){
			if(m_bFree){
				delete [](*m_ppMem);
				m_bFree = 0;
			}
			free(m_ppMem);
			m_ppMem = NULL;
		}
		m_row = 0;
		m_col = 0;
	}

	inline void GetSize(int* pRow, int* pCol){
		*pRow = m_row;
		*pCol = m_col;
	}
	inline T* Get1DPointer(){
		if(NULL==m_ppMem)
			return NULL;
		return *m_ppMem;
	}
	inline T** Get2DPointer(){
		return m_ppMem;
	}
public:
	// product of the two matrices
	void prod(TC2DMem<T> & aMtx, const TC2DMem<T> & bMtx) {
		int i, j, k; 
		(*this).Alloc2DMem(aMtx.m_row,bMtx.m_col);
		for (i=0; i<aMtx.m_row; i++)
			for (j = 0; j<bMtx.m_col; j++){
				(*this).m_ppMem[i][j]=0;
				for (k=0; k<aMtx.m_col; k++)
					(*this).m_ppMem[i][j] += aMtx[i][k] * bMtx[k][j];
			}
		//return (*this);
	}
	// eigen explore of the two matrices
	inline bool eigen2x2(TC2DMem<T> & val, TC2DMem<T> & vec) {
		val.Alloc2DMem(2,2);
		vec.Alloc2DMem(2,2);
		memset(val.Get1DPointer(),0,sizeof(T)*4);
		memset(vec.Get1DPointer(),0,sizeof(T)*4);
		const double a11 = (*this).m_ppMem[0][0], a12 = (*this).m_ppMem[0][1], a21 = (*this).m_ppMem[1][0], a22 = (*this).m_ppMem[1][1];
		const double tr		= a11+a22;
		const double det	= a11*a22 - a12*a21;
		double f = tr*tr - (4. * det);
		if (f<0) return false;
		//if (f<0) cimg::warn("CImg<%s>::eigen() : Complex eigenvalues",pixel_type());
		f = sqrt(f);
		val[0][0]=(T)(tr-f)*.5;
		val[1][1]=(T)(tr+f)*.5;

		f = sqrt(a12*a12 + (val[0][0]-a11)*(val[0][0]-a11) /*+ 0.00000001*/);
		vec[0][0] = (T)   a12/f;
		vec[0][1] = (T) (val[0][0]-a11)/f;
		f = sqrt(a12*a12 + (val[1][1]-a11)*(val[1][1]-a11) /*+ 0.00000001*/);
		vec[1][0] = (T) a12/f;
		vec[1][1] = (T)(val[1][1]-a11)/f;
		return true;
		//return (*this);
	}
	// compute the inverse matrix of the 2x2 matrix
 	inline void inverse2x2(TC2DMem<T> & mat) const{
		mat.Alloc2DMem(2,2);
		memset(mat.Get1DPointer(),0,sizeof(T)*4);
		const double a11 = (*this).m_ppMem[0][0], a12 = (*this).m_ppMem[0][1], a21 = (*this).m_ppMem[1][0], a22 = (*this).m_ppMem[1][1];
		const double inv_det	= 1./(a11*a22 - a12*a21);
		mat[0][0] = inv_det * a22;
		mat[0][1] = inv_det * -1. * a12;
		mat[1][0] = inv_det * -1. * a21;
		mat[1][1] = inv_det * a11;
	}
	//compute the convolution matrix of the 2x2 matrix
	inline void conv2(TC2DMem<T> & horizontalFilter, TC2DMem<T> & verticalFilter, TC2DMem<T> & filteredMatrix) {
		// exception handle
		try{
			if((horizontalFilter.m_col == 0) || (this->m_col <= horizontalFilter.m_col)){
				throw "illegal horizontal filter.";
			}
			if((verticalFilter.m_row == 0) || (this->m_row <= verticalFilter.m_row)){
				throw "illegal vertical filter.";
			}
		}catch (char* e){				
			printf("warning: %s\nPlease try again!",e);
		}		
		
		TC2DMem<T> tempMatrix;
		tempMatrix.Alloc2DMem(this->m_row,this->m_col);
		memset(tempMatrix.Get1DPointer(),0,sizeof(T)*tempMatrix.m_row*tempMatrix.m_col);	
		filteredMatrix.Alloc2DMem(this->m_row, this->m_col);
		memset(filteredMatrix.Get1DPointer(),0,sizeof(T)*filteredMatrix.m_row*filteredMatrix.m_col);		
		
		//horizontal filtering
		for(int row = 0; row < this->m_row; row++){
			for(int column = 0; column < this->m_col-horizontalFilter.m_col+1; column++){				
				for(int ele = 0; ele<horizontalFilter.m_col; ele++){					
					tempMatrix[row][column+horizontalFilter.m_col/2] += (*this)[row][column+ele]*horizontalFilter[0][ele];
				}							
			}
			for(int b = horizontalFilter.m_col/2, c = b+1, a= b-1; b>=1; b--,c--,a--){
				tempMatrix[row][a] = 2*tempMatrix[row][b]-tempMatrix[row][c];
			}
			for(int c = this->m_col-horizontalFilter.m_col/2, b = c-1, a = c-2; c<this->m_col; c++, b++, a++){
				tempMatrix[row][c] = 2*tempMatrix[row][b]-tempMatrix[row][a];
			}
		}	
		
		//vertical filtering
		for(int column = 0; column < this->m_col; column++){
			for(int row = 0; row < this->m_row-verticalFilter.m_row+1; row++){				
				for(int ele = 0; ele<verticalFilter.m_row; ele++){
					 filteredMatrix[row+verticalFilter.m_row/2][column] += tempMatrix[row+ele][column]*verticalFilter[ele][0];
				}				
			}
			for(int b = verticalFilter.m_row/2, c = b+1, a= b-1; b>=1; b--,c--,a--){
				filteredMatrix[a][column] = 2*filteredMatrix[b][column]-filteredMatrix[c][column];
			}
			for(int c = this->m_row-verticalFilter.m_row/2, b = c-1, a = c-2; c<this->m_row; c++, b++, a++){
				filteredMatrix[c][column] = 2*filteredMatrix[b][column]-filteredMatrix[a][column];
			}
		}
		
	}


public:
	inline T* operator[](int nRow) const{
		return m_ppMem[nRow];
	}

	//inline T* operator*(T* mat) const{	
	//	return *this;
	//}
	inline void operator*=(const float d) {
		for (int i=0;i<m_row;i++)
			for (int j = 0; j<m_col; j++)
				this->m_ppMem[i][j] *= d; 
		return;
	}

	inline void operator^(const float power) {
		for (int i=0;i<m_row;i++)
			for (int j = 0; j<m_col; j++)
				this->m_ppMem[i][j] = pow(this->m_ppMem[i][j], power);
				//sqrt(this->m_ppMem[i][j]); 
		return;
	}
};
/*
template<class T>
TC2DMem<T>::TC2DMem(){
	m_row = 0;
	m_col = 0;
	m_ppMem = NULL;
}


template<class T>	
TC2DMem<T>::~TC2DMem(){
	Free2DMem();
}

template<class T>
T** TC2DMem<T>::Alloc2DMem(int nRow, int nCol){
	T** ppMem;
	T* pMem;

	T** ppCur;
	T* pCur;

	int tmp;

	if(nRow<1)
		return NULL;
	if(nCol<1)
		return NULL;
	m_row = nRow;
	m_col = nCol;
	Free2DMem();

	ppMem = (T**)malloc(m_row*sizeof(T*));
	pMem = new T[m_row*m_col];
	ppCur = ppMem;
	pCur = pMem;

	for(tmp=0; tmp<m_row; tmp++){
		*ppCur = pCur;
		pCur += m_col;
		ppCur++;
	}

	m_ppMem = ppMem;

	return m_ppMem;
}

template<class T>
void TC2DMem<T>::Free2DMem(){
	if(NULL!=(m_ppMem)){
		delete [](*m_ppMem);
		free(m_ppMem);
		m_ppMem = NULL;
	}
}

template<class T>
void TC2DMem<T>::GetSize(int* pRow, int* pCol){
	*pRow = m_row;
	*pCol = m_col;
}

template<class T>
T* TC2DMem<T>::Get1DPointer(){
	if(NULL==m_ppMem)
		return NULL;
	return *m_ppMem;
}

template<class T>
T** TC2DMem<T>::Get2DPointer(){
	return m_ppMem;
}
*/
#endif