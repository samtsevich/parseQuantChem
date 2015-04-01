/********************************************************************************
����� CMatrix, ��������������� �������� �������� � ���������.
������ ���� ���������������� �� �������� AS IS.
����������� ������������� ����������� ��� ������� ���������� ������� ���������.

�����������: ������ aka CyBOSSeR �����
E-mail: cybosser@gmail.com
********************************************************************************/

#ifndef _MATRIX_
#define _MATRIX_

#include <memory.h>

class CMatrix
{
private:
	class CRow
	{
	private:
		friend class CMatrix;

		CRow(int* pRow, const int nColCount);
		// ��������� �� ������ ������ � ������
		int* m_pFirstCellInRow;
		// ���������� ����� � ������
		int m_nRowElementsCount;
	public:
		int &operator[](const int nCellInRow) const;
	};
public:
	CMatrix(const int nRowCount = 0, const int nColCount = 0);
	CMatrix(const CMatrix &aMatrix);
	~CMatrix(void);

	CRow operator[](const int nRow) const;
	CMatrix operator*(const CMatrix &aMatrix) const;
	CMatrix operator*(const int n) const;
	CMatrix operator+(const CMatrix &aMatrix) const;
	CMatrix operator-(const CMatrix &aMatrix) const;
	CMatrix operator^(const int n) const;
	
	CMatrix& operator=(const CMatrix &aMatrix);
	CMatrix& operator*=(const int n);
	CMatrix& operator*=(const CMatrix &aMatrix);
	CMatrix& operator+=(const CMatrix &aMatrix);
	CMatrix& operator-=(const CMatrix &aMatrix);
	CMatrix& operator^=(const int n);

	int GetRowCount(void) const;
	int GetColCount(void) const;
	void Resize(const int nNewRowCount, const int nNewColCount);
private:
	// ���������� ����� � �������
	int m_nRowCount;
	// ���������� �������� � �������
	int m_nColCount;
	// ������ ���������� �������� �������
	int *m_pArray;
};
//===============================================================================
// �������-����� ������ CMatrix.
//===============================================================================
CMatrix::CMatrix(const int nRowCount/*=0*/, const int nColCount/*=0*/)
: m_nRowCount(nRowCount)
, m_nColCount(nColCount)
, m_pArray(new int[nRowCount*nColCount])
{ memset(m_pArray, 0, sizeof(int)*m_nRowCount*m_nColCount); }
//-------------------------------------------------------------------------------
CMatrix::CMatrix(const CMatrix &aMatrix)
: m_nRowCount(aMatrix.m_nRowCount)
, m_nColCount(aMatrix.m_nColCount)
, m_pArray(new int[aMatrix.m_nRowCount*aMatrix.m_nColCount])
{ memcpy(m_pArray, aMatrix.m_pArray, sizeof(int)*m_nRowCount*m_nColCount);}
//-------------------------------------------------------------------------------
CMatrix::~CMatrix(void)
{
	if(m_pArray)
		delete [] m_pArray;
}
//-------------------------------------------------------------------------------
// �������-����, ������������ ������ ����������� ����� CRow, ��������� ���
// ������� � ������������ ������ ��������� ������� nRow.
//-------------------------------------------------------------------------------
CMatrix::CRow CMatrix::operator[](const int nRow) const
{ return CRow(m_pArray+nRow*m_nColCount, m_nColCount); }
//-------------------------------------------------------------------------------
// �������-���� ������������ ���� ������.
// � ������, ���� ���������� �������� ������ ������� �� ����� ����� ����� ������
// �������, ������������ ������ ������� �������� 0x0.
//-------------------------------------------------------------------------------
CMatrix CMatrix::operator*(const CMatrix &aMatrix) const
{
	CMatrix tmpMatrix(*this);

	tmpMatrix *= aMatrix;

	return tmpMatrix;
}
//-------------------------------------------------------------------------------
// �������-���� ��������� ������� �� �������� ����� n.
//-------------------------------------------------------------------------------
CMatrix CMatrix::operator*(const int n) const
{
	CMatrix tmpMatrix(*this);

	tmpMatrix *= n;

	return tmpMatrix;
}
//-------------------------------------------------------------------------------
// �������-���� �������� ���� ������.
//-------------------------------------------------------------------------------
CMatrix CMatrix::operator+(const CMatrix &aMatrix) const
{
	CMatrix tmpMatrix(*this);

	tmpMatrix += aMatrix;

	return tmpMatrix;
}
//-------------------------------------------------------------------------------
// �������-���� �������� ���� ������.
//-------------------------------------------------------------------------------
CMatrix CMatrix::operator-(const CMatrix &aMatrix) const
{
	CMatrix tmpMatrix(*this);

	tmpMatrix -= aMatrix;

	return tmpMatrix;
}
//-------------------------------------------------------------------------------
// �������-���� ���������� ������� � ����������� ������� n.
//-------------------------------------------------------------------------------
CMatrix CMatrix::operator^(const int n) const
{
	CMatrix tmpMatrix(*this);

	tmpMatrix ^= n;

	return tmpMatrix;
}
//-------------------------------------------------------------------------------
// �������-���� ���������� ������� �� ������ �������.
// � ������, ���� ���������� �������� ������ ������� �� ����� ����� ����� ������
// �������, �������� ������� ������������� � ������� �������� 0x0.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator*=(const CMatrix &aMatrix)
{
	if(m_nColCount != aMatrix.m_nRowCount)
	{
		Resize(0, 0);
		return *this;
	}

	CMatrix tmpMatrix(m_nRowCount, aMatrix.m_nColCount);

	for(int i=0; i<tmpMatrix.m_nRowCount; i++)
		for(int j=0; j<tmpMatrix.m_nColCount; j++)
			for(int k=0; k<m_nColCount; k++)
				tmpMatrix[i][j]+=(*this)[i][k]*aMatrix[k][j];

	(*this) = tmpMatrix;

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� ���������� ������� �� �������� ����� n.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator*=(const int n)
{
	int nLenghtOfArray = m_nRowCount * m_nColCount;

	for(int i=0; i<nLenghtOfArray; i++)
		m_pArray[i] *= n;

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� �������� ������� � ������ ��������.
// � ������, ���� ������� ����� ��������� �����������, ������� �����������
// �������� ������� � ������� ������������ 0x0.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator+=(const CMatrix &aMatrix)
{
	if(m_nRowCount != aMatrix.m_nRowCount || m_nColCount != aMatrix.m_nColCount)
	{
		Resize(0, 0);
		return *this;
	}

	int nLenghtOfArray = m_nRowCount * m_nColCount;

	for(int i=0; i<nLenghtOfArray; i++)
		m_pArray[i] += aMatrix.m_pArray[i];

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� ��������� ������� aMatrix �� �������.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator-=(const CMatrix &aMatrix)
{
	if(m_nRowCount != aMatrix.m_nRowCount || m_nColCount != aMatrix.m_nColCount)
	{
		Resize(0,0);
		return *this;
	}

	int nLenghtOfArray = m_nRowCount * m_nColCount;

	for(int i=0; i<nLenghtOfArray; i++)
		m_pArray[i] -= aMatrix.m_pArray[i];

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� ���������� ������� � ����������� ������� n.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator^=(const int n)
{
	CMatrix tmpMatrix(*this);

	for(int i=1; i<n; i++)
		(*this) *= tmpMatrix;

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� ������������.
//-------------------------------------------------------------------------------
CMatrix& CMatrix::operator=(const CMatrix &aMatrix)
{
	if(this==&aMatrix)
		return *this;

	if(m_pArray)
		delete [] m_pArray;

	m_nRowCount = aMatrix.m_nRowCount;
	m_nColCount = aMatrix.m_nColCount;

	m_pArray = new int[m_nRowCount*m_nColCount];

	memcpy(m_pArray, aMatrix.m_pArray, sizeof(int)*m_nRowCount*m_nColCount);

	return *this;
}
//-------------------------------------------------------------------------------
// �������-���� ��������� ���������� ����� �������.
//-------------------------------------------------------------------------------
inline int CMatrix::GetRowCount(void) const
{ return m_nRowCount; }
//-------------------------------------------------------------------------------
// �������-���� ��������� ���������� �������� �������.
//-------------------------------------------------------------------------------
inline int CMatrix::GetColCount(void) const
{ return m_nColCount; }
//-------------------------------------------------------------------------------
// �������-���� ��������� �������� �������.
//-------------------------------------------------------------------------------
void CMatrix::Resize(const int nNewRowCount, const int nNewColCount)
{
	if(m_pArray)
		delete [] m_pArray;

	m_nRowCount = nNewRowCount;
	m_nColCount = nNewColCount;

	m_pArray = new int[m_nRowCount*m_nColCount];
	
	memset(m_pArray, 0, sizeof(int)*m_nRowCount*m_nColCount);
}
//===============================================================================
// �������-����� ����� CMatrix::CRow.
//===============================================================================
CMatrix::CRow::CRow(int *pRow, const int nColCount)
: m_pFirstCellInRow(pRow)
, m_nRowElementsCount(nColCount)
{}
//-------------------------------------------------------------------------------
// �������� ������������ ��������� ������ � ������. 
//-------------------------------------------------------------------------------
int &CMatrix::CRow::operator [](const int nCellInRow) const
{ 
	return (nCellInRow < m_nRowElementsCount)? m_pFirstCellInRow[nCellInRow]:
		                                       m_pFirstCellInRow[0];
}
//-------------------------------------------------------------------------------

#endif