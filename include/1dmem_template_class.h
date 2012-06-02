#ifndef __1DMEM_TEMPLATE_CLASS_H__
#define __1DMEM_TEMPLATE_CLASS_H__

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

#include <stdlib.h>

template<class T>
class TC1DMem{
public:
	int m_size;
	T* m_pMem;
private:
	bool m_bFree;
public:
	TC1DMem(){
		m_size = 0;
		m_bFree = 0;
		m_pMem = NULL;
	}
	~TC1DMem(){
		Free1DMem();
	}
	T* Alloc1DMem(int nSize){
		Free1DMem();
		if(nSize<1)
			return NULL;
		m_size = nSize;

		m_pMem = new T[m_size];
		m_bFree = 1;
		return m_pMem;
	}
	T* AssignMemory(int nSize, T* pMem){
		Free1DMem();
		if(nSize<1)
			return NULL;
		m_size = nSize;

		m_pMem = pMem;
		return m_pMem;
	}
	void Free1DMem(){
		if(NULL!=(m_pMem)){
			if(m_bFree){
				delete [](m_pMem);
				m_bFree = 0;
			}
			m_pMem = NULL;
		}
		m_size = 0;
	}

	inline void GetSize(int* pSize){
		*pSize = m_size;
	}
	inline T* GetPointer(){
		return m_pMem;
	}
public:
	inline T& operator[](int nElement) const{
		return m_pMem[nElement];
	}
};

#endif